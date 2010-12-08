#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_H

/// Crystallographic least-squares

#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/math/accumulators.h>
#include <scitbx/sparse/matrix.h>

#include <cctbx/sgtbx/seminvariant.h>

#include <smtbx/refinement/constraints/scatterer_parameters.h>
#include <smtbx/structure_factors/direct/standard_xray.h>

#include <algorithm>
#include <iterator>


namespace smtbx { namespace refinement { namespace least_squares {

  namespace lstbx = scitbx::lstbx;

  /** \brief Restraints to prevent the structure to freely move as a whole
   along the continuous shift directions if any.

   The implementation has a key limitation: the continuous shift directions
   must lie along unit cell axes if there are atoms on special positions on
   axes parallel to continuous shift directions.
   It will therefore work for all space-groups in standard settings
   but it will not work in, e.g., a setting of R3 with
   the 3-fold axis along the direction (1,1,1), as in the space group
   Hall: R 3 (-y+z, x+z, -x+y+z), with an atom at the position (0.5, 0.5, 0.5).
   As far as we know, this is a limitation of the programs ShelXL
   and Crystals too. In order to keep it simple, an exception is thrown if
   a continuous shift directions is found to be non-trivial, whether there are
   atoms on the wrong kind of special positions or not.

   The weight of the restraints is dynamically adjusted to the normal matrix
   they are added to, as a relative weight times the smallest diagonal element.
   */
  template <typename FloatType>
  class floating_origin_restraints
  {
  public:
    typedef FloatType scalar_t;

    scalar_t floating_origin_restraint_relative_weight;

    af::small<af::shared<scalar_t>, 3> singular_directions;
    af::small<scitbx::vec3<scalar_t>, 3> origin_shifts;

  public:
    /// Construct with the given relative weight
    /** It precomputes the singular vectors of the normal matrix
        for the independent parameters after performing the reparametrisations
        dictated by the constraints. The necessary information about the latter
        is passed in the transpose of the Jacobian matrix.
     */
    floating_origin_restraints(
      sgtbx::space_group const &space_group,
      af::const_ref<constraints::scatterer_parameters> const &params,
      scitbx::sparse::matrix<FloatType> const
      &jacobian_transpose_matching_grad_fc,
      scalar_t floating_origin_restraint_relative_weight)

      : floating_origin_restraint_relative_weight(
          floating_origin_restraint_relative_weight)
    {
      // Floating origin restraints: pre-compute singular directions
      sgtbx::structure_seminvariants seminvariants(space_group);
      SMTBX_ASSERT(seminvariants.continuous_shifts_are_principal());
      af::small<sgtbx::ss_vec_mod, 3> const
      &vm = seminvariants.vectors_and_moduli();
      for (int i=0; i<vm.size(); ++i) {
        if (vm[i].m != 0) continue;
        // allowed continuous origin shift:
        scitbx::vec3<scalar_t> v(vm[i].v[0], vm[i].v[1], vm[i].v[2]);
        origin_shifts.push_back(v);

        // Fill singular_dir with the singular direction corresponding to v
        af::shared<scalar_t> singular_dir(af::reserve(5*params.size()));
        std::back_insert_iterator< af::shared<scalar_t> > g(singular_dir);
        for (std::size_t i=0; i<params.size(); ++i) {
          BOOST_FOREACH (constraints::asu_parameter const *p,
                         params[i].ordered())
          {
            if (p == params[i].site) std::copy(v.begin(), v.end(), g);
            else std::fill_n(g, p->size(), scalar_t(0));
          }
        }

        // Save it
        singular_directions.push_back(
          jacobian_transpose_matching_grad_fc*singular_dir.ref());
      }
    }

    /// Add floating origin restraints to the given normal equations
    void add_to(lstbx::normal_equations::linear_ls<scalar_t> &normal_eqns) {
      if (!floating_origin_restraint_relative_weight) return;
      af::ref<scalar_t, af::packed_u_accessor>
      a = normal_eqns.normal_matrix().ref();
      scitbx::math::accumulator::min_max_accumulator<scalar_t> acc(a(0,0));
      for (int i=1; i<a.n_rows(); ++i) acc(a(i,i));
      scalar_t w = floating_origin_restraint_relative_weight * acc.max();
      for (int i=0; i<singular_directions.size(); ++i) {
        normal_eqns.add_equation(0, singular_directions[i].const_ref(), w);
      }
    }
  };


  /** \brief Build normal equations for the given data, model, weighting
       and constraints.

      The constraints is performed with a reparametrisation whose Jacobian
      transpose is passed as an argument.
   */
  template <typename FloatType>
  struct build_normal_equations
  {
    //! Default constructor. Some data members are not initialized!
    build_normal_equations() {}

    template <template<typename> class NormalEquations,
              template<typename> class WeightingScheme,
              class OneMillerIndexFcalc>
    build_normal_equations(
      NormalEquations<FloatType> &normal_equations,
      af::const_ref<miller::index<> > const &miller_indices,
      af::const_ref<FloatType> const &data,
      af::const_ref<FloatType> const &sigmas,
      af::const_ref<std::complex<FloatType> > const &f_mask,
      WeightingScheme<FloatType> const &weighting_scheme,
      FloatType scale_factor,
      OneMillerIndexFcalc &f_calc_function,
      scitbx::sparse::matrix<FloatType> const
        &jacobian_transpose_matching_grad_fc,
      bool objective_only=false)
    :
      f_calc_(miller_indices.size()),
      weights_(miller_indices.size())
    {
      // Accumulate equations Fo(h) ~ Fc(h)
      SMTBX_ASSERT(miller_indices.size() == data.size())
                  (miller_indices.size())(data.size());
      SMTBX_ASSERT(data.size() == sigmas.size())(data.size())(sigmas.size());
      if (f_mask.size()){
        SMTBX_ASSERT(f_mask.size() == f_mask.size())(data.size())(data.size());
      }
      if (objective_only) {
        for (int i_h=0; i_h<miller_indices.size(); ++i_h) {
          miller::index<> const &h = miller_indices[i_h];
          if (f_mask.size()) {
            f_calc_function.evaluate(h, f_mask[i_h]);
          }
          else {
            f_calc_function.evaluate(h);
          }
          FloatType observable = f_calc_function.observable;
          FloatType weight = weighting_scheme(data[i_h], sigmas[i_h],
                                              observable, scale_factor);
          f_calc_[i_h] = f_calc_function.f_calc;
          weights_[i_h] = weight;
          normal_equations.add_residual(observable,
                                        data[i_h],
                                        weight);
        }
        normal_equations.finalise(/*objective_only=*/true);
      }
      else {
        for (int i_h=0; i_h<miller_indices.size(); ++i_h) {
          miller::index<> const &h = miller_indices[i_h];
          if (f_mask.size()) {
            f_calc_function.linearise(h, f_mask[i_h]);
          }
          else {
            f_calc_function.linearise(h);
          }
          FloatType observable = f_calc_function.observable;
          af::shared<FloatType> gradient =
            jacobian_transpose_matching_grad_fc*f_calc_function.grad_observable;
          FloatType weight = weighting_scheme(data[i_h], sigmas[i_h],
                                              observable, scale_factor);
          f_calc_[i_h] = f_calc_function.f_calc;
          weights_[i_h] = weight;
          normal_equations.add_equation(observable,
                                        gradient.ref(),
                                        data[i_h],
                                        weight);
        }
        normal_equations.finalise();
      }
    }

    af::shared<std::complex<FloatType> > f_calc() { return f_calc_; }

    af::shared<FloatType> weights() { return weights_; }

  private:
    af::shared<std::complex<FloatType> > f_calc_;
    af::shared<FloatType> weights_;

  };

}}}


#endif // GUARD
