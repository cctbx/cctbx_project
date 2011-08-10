#ifndef SMTBX_REFINEMENT_RESTRAINTS
#define SMTBX_REFINEMENT_RESTRAINTS

#include <scitbx/math/accumulators.h>
#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>

#include <cctbx/sgtbx/seminvariant.h>

#include <smtbx/refinement/constraints/scatterer_parameters.h>

#include <iterator>


namespace smtbx { namespace refinement { namespace restraints {

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

   If \f$A\f$ is the original normal matrix, the restraint replace \f$A\f$
   by \f$A + \mu I\f$. The weight \f$\mu\f$ of the restraints is dynamically
   adjusted to \f$A\f$, so as to fullfill two requirements:

   (a) \f$A + \mu I\f$ is non-singular enough that the L.S. problem can
   be reliably solved; and

   (b) \f$(A + \mu I)^{-1}\f$ is approximately the pseudo-inverse of \f$A\f$.

   Requirement (b) is essential to obtaining faithful parameter esd's and
   correlations.
   */
  template <typename FloatType>
  class origin_fixing
  {
  public:
    typedef FloatType scalar_t;

    scalar_t relative_weight;

    /// Singular directions in the space of crystallographic parameters
    af::small<af::shared<scalar_t>, 3> singular_directions;

    af::small<scitbx::vec3<scalar_t>, 3> origin_shifts;

    /// A mask pinpointing which parameters span the singular space
    af::shared<bool> singular_space;

  public:
    /// Construct with the given relative weight
    /** It precomputes the singular vectors of the normal matrix
        for the crystallographic parameters.
     */
    origin_fixing(
      sgtbx::space_group const &space_group,
      af::const_ref<constraints::scatterer_parameters> const &params,
      scalar_t relative_weight)

      : relative_weight(relative_weight)
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
        singular_directions.push_back(singular_dir);
      }

      // Pinpoint coordinates participating in the singularity
      if (singular_directions.size()) {
        int n = singular_directions[0].size();
        singular_space = af::shared<bool>(n, false);
        for (int i=0; i<n; ++i) {
          for (int j=0; j<singular_directions.size(); ++j) {
            singular_space[i] = singular_space[i] || singular_directions[j][i];
          }
        }
      }
    }

    /// Add floating origin restraints to the given normal equations
    /** The Jacobian is that of the reparametrisation of crystallographic
        parameters as function of the independent parameters that the
        normal equations are built with
     */
    void add_to(lstbx::normal_equations::linear_ls<scalar_t> &normal_eqns,
                scitbx::sparse::matrix<FloatType> const
                &jacobian_transpose_matching_grad_fc)
    {
      if (!relative_weight) return;
      if (!singular_space.size()) return;
      af::ref<scalar_t, af::packed_u_accessor>
      a = normal_eqns.normal_matrix().ref();
      scitbx::math::accumulator::min_max_accumulator<scalar_t> acc(a(0,0));
      for (int i=1; i<a.n_rows(); ++i) {
        if (singular_space[i]) acc(a(i,i));
      }
      scalar_t w = relative_weight * acc.max();
      for (int i=0; i<singular_directions.size(); ++i) {
        af::shared<scalar_t> reparametrised_singular_dir =
        jacobian_transpose_matching_grad_fc*singular_directions[i].const_ref();
        normal_eqns.add_equation(0, reparametrised_singular_dir.const_ref(), w);
      }
    }
  };

}}}



#endif // Guard
