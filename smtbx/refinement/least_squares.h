#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/shared_reductions.h>
#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/matrix/cholesky.h>
#include <scitbx/math/accumulators.h>

#include <cctbx/math/cos_sin_table.h>
#include <cctbx/xray/scatterer.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/coordinates.h>
#include <cctbx/adptbx.h>

#include <smtbx/error.h>
#include <smtbx/structure_factors/direct/standard_xray.h>

#include <algorithm>
#include <iterator>


namespace smtbx { namespace refinement { namespace least_squares {

  /** \brief Restraints to prevent the structure to freely move as a whole
   along the continuous shift directions if any.

   The weight of the restraints is dynamically adjusted to the normal matrix
   they are added to, as a relative weight times the largest diagonal element.

   The implementation is sub-optimal but computations are inexpensive
   compared to the rest of a L.S. refinement cycle. Thus, not worth optimising.

   */
  template <typename FloatType>
  class floating_origin_restraints
  {
  public:
    typedef FloatType scalar_t;
    typedef xray::scatterer<scalar_t> scatterer_t;

    scalar_t floating_origin_restraint_relative_weight;

    af::small<af::shared<scalar_t>, 3> singular_directions;
    af::small<scitbx::vec3<scalar_t>, 3> origin_shifts;

  public:
    /// Construct with the given relative weight
    /** The singular vectors of the crystallographic normal matrix are
     precomputed, taking into account special position constraints.
     */
    floating_origin_restraints(
      sgtbx::space_group const &space_group,
      sgtbx::site_symmetry_table const &site_symmetry_table,
      af::shared<scatterer_t> const &scatterers,
      scalar_t floating_origin_restraint_relative_weight)

      : floating_origin_restraint_relative_weight(
          floating_origin_restraint_relative_weight)
    {
      // Floating origin restraints: pre-compute singular directions
      sgtbx::structure_seminvariants seminvariants(space_group);
      af::small<sgtbx::ss_vec_mod, 3> const
      &vm = seminvariants.vectors_and_moduli();
      for (int i=0; i<vm.size(); ++i) {
        if (vm[i].m != 0) continue;
        // allowed continuous origin shift:
        scitbx::vec3<scalar_t> v(vm[i].v[0], vm[i].v[1], vm[i].v[2]);
        origin_shifts.push_back(v);

        // Fill singular_dir with the singular direction corresponding to v
        af::shared<scalar_t> singular_dir(af::reserve(3*scatterers.size()));
        std::back_insert_iterator< af::shared<scalar_t> > g(singular_dir);
        for (int i_sc=0; i_sc<scatterers.size(); ++i_sc) {
          scatterer_t const& sc = scatterers[i_sc];
          const sgtbx::site_symmetry_ops& op = site_symmetry_table.get(i_sc);
          if (sc.flags.grad_site()) {
            if (op.is_point_group_1()) {
              std::copy(v.begin(), v.begin() + 3, g);
            }
            else {
              sgtbx::site_constraints<scalar_t> const
              &site_c = op.site_constraints();
              af::small<scalar_t,3> v1 = site_c.independent_params(v);
              std::copy(v1.begin(), v1.end(), g);
            }
          }
          if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) {
            *g = 0;
          }
          if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
            std::fill_n(g, op.adp_constraints().n_independent_params(),
                        scalar_t(0));
          }
          if (sc.flags.grad_occupancy()) {
            *g = 0;
          }
        }

        // Save it
        singular_directions.push_back(singular_dir);
      }
    }

    /// Add floating origin restraints to the given normal equations
    void add_to(scitbx::lstbx::normal_equations<scalar_t> &normal_eqns) {
      if (!floating_origin_restraint_relative_weight) return;
      typedef scitbx::lstbx::normal_equations<scalar_t> normal_eqns_t;
      af::ref<scalar_t, af::packed_u_accessor>
      a = normal_eqns.normal_matrix().ref();
      scitbx::math::accumulator::min_max_accumulator<scalar_t> acc(a(0,0));
      for (int i=1; i<a.n_rows(); ++i) acc(a(i,i));
      scalar_t w = floating_origin_restraint_relative_weight * acc.max();
      for (int i=0; i<singular_directions.size(); ++i) {
        normal_eqns.add_equation(0,
                                 singular_directions[i].const_ref(),
                                 w);
      }
    }
  };


  /** \brief Build normal equations for the given data, model and weighting
   enforcing the given constraints on refined parameters.
   */
  template <typename FloatType,
            template<typename> class NormalEquations,
            template<typename> class WeightingScheme,
            class OneMillerIndexLinearisation,
            class ConstraintsType>
  void
  build_normal_equations(NormalEquations<FloatType> &normal_equations,
                         af::const_ref<miller::index<> > const &miller_indices,
                         af::const_ref<FloatType> const &data,
                         af::const_ref<FloatType> const &sigmas,
                         af::const_ref<std::complex<FloatType> > const &f_mask,
                         WeightingScheme<FloatType> const &weighting_scheme,
                         FloatType scale_factor,
                         OneMillerIndexLinearisation &one_h_linearisation,
                         ConstraintsType &constraints)
  {
    // Accumulate equations Fo(h) ~ Fc(h)
    SMTBX_ASSERT(miller_indices.size() == data.size())
                (miller_indices.size())(data.size());
    SMTBX_ASSERT(data.size() == sigmas.size())(data.size())(sigmas.size());
    if (f_mask.size()){
      SMTBX_ASSERT(f_mask.size() == f_mask.size())(data.size())(data.size());
    }
    for (int i_h=0; i_h<miller_indices.size(); ++i_h) {
      miller::index<> const &h = miller_indices[i_h];
      if (f_mask.size()) {
        one_h_linearisation.compute(h, f_mask[i_h]);
      }
      else {
        one_h_linearisation.compute(h);
      }
      constraints.apply_chain_rule(one_h_linearisation.grad_observable);
      FloatType observable = one_h_linearisation.observable;
      FloatType weight = weighting_scheme(data[i_h], sigmas[i_h],
                                          scale_factor * observable);
      normal_equations.add_equation(observable,
                                    constraints.independent_gradients().begin(),
                                    data[i_h],
                                    weight);
    }
    normal_equations.finalise();
  }

  /// Special position constraints
  template <typename FloatType>
  class special_position_constraints
  {
  public:
    typedef FloatType scalar_t;
    typedef xray::scatterer<scalar_t> scatterer_t;

    std::size_t n_independent_params, n_crystallographic_params;

  private:
    uctbx::unit_cell const &unit_cell;
    sgtbx::site_symmetry_table const &site_symmetry_table;
    af::ref_owning_shared<scatterer_t> scatterers;
    af::shared<scalar_t> independent_grads;

  public:
    special_position_constraints(uctbx::unit_cell const &unit_cell,
                                 sgtbx::site_symmetry_table const
                                 &site_symmetry_table,
                                 af::shared<scatterer_t> const &scatterers)
      : unit_cell(unit_cell),
        site_symmetry_table(site_symmetry_table),
        scatterers(scatterers),
        n_independent_params(0), n_crystallographic_params(0)
    {
      for (int i_sc=0; i_sc<scatterers.size(); ++i_sc) {
        scatterer_t const& sc = scatterers[i_sc];
        const sgtbx::site_symmetry_ops& op = site_symmetry_table.get(i_sc);
        if (sc.flags.grad_site()) {
          n_crystallographic_params += 3;
          n_independent_params += op.site_constraints().n_independent_params();
        }
        if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) {
          n_crystallographic_params += 1;
          n_independent_params += 1;
        }
        if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
          n_crystallographic_params += 6;
          n_independent_params += op.adp_constraints().n_independent_params();
        }
        if (sc.flags.grad_occupancy()) {
          n_crystallographic_params += 1;
          n_independent_params += 1;
        }
      }
      independent_grads.resize(n_independent_params);
    }

    af::shared<scalar_t> const independent_gradients() const {
      return independent_grads;
    }

    void
    apply_chain_rule(af::const_ref<scalar_t> const &crystallographic_gradients)
    {
      if (!independent_grads.size()) return;
      scalar_t const *xg = crystallographic_gradients.begin();
      scalar_t       *g  = independent_grads.begin();
      for (int i_sc=0; i_sc<scatterers.size(); ++i_sc) {
        scatterer_t const& sc = scatterers[i_sc];
        const sgtbx::site_symmetry_ops& op = site_symmetry_table.get(i_sc);
        if (sc.flags.grad_site()) {
          if (op.is_point_group_1()) {
            g = std::copy(xg, xg+3, g);
          }
          else {
            scitbx::vec3<scalar_t> grad_site(xg);
            sgtbx::site_constraints<scalar_t> const
            &site_c = op.site_constraints();
            af::small<scalar_t,3>
            ind_grad_site = site_c.independent_gradients(grad_site);
            g = std::copy(ind_grad_site.begin(), ind_grad_site.end(), g);
          }
          xg += 3;
        }
        if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) {
          *g++ = *xg++;
        }
        if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
          scitbx::sym_mat3<scalar_t>
          frac_grad(xg),
          cart_grad = adptbx::grad_u_star_as_u_cart(unit_cell, frac_grad);
          if (op.is_point_group_1()) {
            g = std::copy(cart_grad.begin(), cart_grad.end(), g);
          }
          else {
            const sgtbx::tensor_rank_2::cartesian_constraints<scalar_t>
            &adp_c = op.cartesian_adp_constraints(unit_cell);
            af::small<scalar_t, 6>
            ind_cart_grad = adp_c.independent_gradients(cart_grad);
            g = std::copy(ind_cart_grad.begin(), ind_cart_grad.end(), g);
          }
          xg += 6;
        }
        if (sc.flags.grad_occupancy()) {
          *g++ = *xg++;
        }
      }
      SCITBX_ASSERT(xg - crystallographic_gradients.begin()
                    == n_crystallographic_params);
      SCITBX_ASSERT(g - independent_grads.begin()
                    == n_independent_params);
    }

    void apply_shifts(af::const_ref<FloatType> const &shifts) const
    {
      scalar_t const *g = shifts.begin();
      for (int i_sc=0; i_sc < scatterers.size(); ++i_sc) {
        scatterer_t &sc = scatterers[i_sc];
        sgtbx::site_symmetry_ops const &op = site_symmetry_table.get(i_sc);
        if (sc.flags.grad_site()) {
          if (op.is_point_group_1()) {
            fractional<scalar_t> delta_site(g);
            sc.site += delta_site;
            g += 3;
          }
          else {
            sgtbx::site_constraints<FloatType> const
            &site_c = op.site_constraints();
            int n = site_c.n_independent_params();
            af::small<scalar_t, 3> shift(g, g+n);
            sc.site += site_c.all_shifts(shift);
            g += n;
          }
        }
        if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) {
          sc.u_iso += *g++;
        }
        if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
          if (op.is_point_group_1()) {
            scitbx::sym_mat3<scalar_t>
            delta_u_cart(g),
            delta_u_star = adptbx::u_cart_as_u_star(unit_cell, delta_u_cart);
            sc.u_star += delta_u_star;
            g += 6;
          }
          else {
            const sgtbx::tensor_rank_2::cartesian_constraints<scalar_t>
            &adp_c = op.cartesian_adp_constraints(unit_cell);
            int n = adp_c.n_independent_parameters();
            af::small<scalar_t, 6> shift(g, g+n);
            scitbx::sym_mat3<scalar_t> delta_u_cart = adp_c.all_params(shift);
            scitbx::sym_mat3<scalar_t>
            delta_u_star = adptbx::u_cart_as_u_star(unit_cell, delta_u_cart);
            sc.u_star += delta_u_star;
            g += n;
          }
        }
        if (sc.flags.grad_occupancy()) {
          sc.occupancy += *g++;
        }
      }
      SCITBX_ASSERT(g - shifts.begin() == shifts.size())
                   (g - shifts.begin())(shifts.size());
    }
  };


}}}


#endif // GUARD
