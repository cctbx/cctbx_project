#ifndef SMTBX_REFINEMENT_RESTRAINTS
#define SMTBX_REFINEMENT_RESTRAINTS

#include <scitbx/vec3.h>
#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>

#include <cctbx/sgtbx/seminvariant.h>

#include <smtbx/refinement/constraints/scatterer_parameters.h>

#include <iterator>


namespace smtbx { namespace refinement { namespace restraints {

  namespace lstbx = scitbx::lstbx;

  /** \brief Restraints to prevent the structure to freely move as a whole
   along the continuous shift directions if any.

   The code uses as much as possible the notations of the reference paper.

   Reference:
   Flack, H. D., and Schwarzenbach, D.
   On the use of least-square restraints for origin fixing in polar space groups.
   Acta Cryst. A 44 (1988), 499Ð506.
   */
  template <typename FloatType>
  class origin_fixing
  {
  public:
    typedef FloatType scalar_t;
    typedef scitbx::vec3<scalar_t> vec3_t;
    typedef constraints::scatterer_parameters::scatterer_type
            scatterer_type;

  private:
    af::small<vec3_t, 3> deltas;
    af::small<af::shared<scalar_t>, 3> whole_deltas;

  public:

    /// Independent directions in which the whole structure may freely float
    af::small<vec3_t, 3> const &origin_shifts() const {
      return deltas;
    }

    /// Whether there are directions in which a structure may freely float
    bool has_floating_directions() const {
      return deltas.size();
    }

    /// Singular direction of the normal matrix for crystallographic parameters
    /** That is before the reparametrisation triggered by constraints
        is applied.
     */
    af::small<af::shared<scalar_t>, 3> const &singular_directions() const {
      return whole_deltas;
    }

    /** \brief Precomputes the origin shifts for the given space group.
     */
    origin_fixing(sgtbx::space_group const &space_group)
    {
      sgtbx::structure_seminvariants seminvariants(space_group);
      af::small<sgtbx::ss_vec_mod, 3> const
      &vm = seminvariants.vectors_and_moduli();
      for (int i=0; i<vm.size(); ++i) {
        if (vm[i].m != 0) continue;
        // allowed continuous origin shift:
        scitbx::vec3<scalar_t> v(vm[i].v[0], vm[i].v[1], vm[i].v[2]);
        deltas.push_back(v);
      }
    }

    ~origin_fixing()
    {}

    /// Scatterer weights
    /** Scatterer weights are the \f$a_m\f$ in the reference paper.
        The i-th element of scatterer_weights shall be the weight
        for the scatterer corresponding to the i-th element of params.
        This member function will be called from member function add_to.
     */
    virtual af::shared<scalar_t>
    weights(lstbx::normal_equations::linear_ls<scalar_t> &normal_eqns,
            scitbx::sparse::matrix<scalar_t> const
            &jacobian_transpose_matching_grad_fc,
            af::shared<constraints::scatterer_parameters> const &params) = 0;

    /// Add floating origin restraints to the given normal equations
    /** The Jacobian is that of the reparametrisation of crystallographic
        parameters as function of the independent parameters that the
        normal equations are built with.
     */
    void add_to(lstbx::normal_equations::linear_ls<scalar_t> &normal_eqns,
                scitbx::sparse::matrix<scalar_t> const
                &jacobian_transpose_matching_grad_fc,
                af::shared<constraints::scatterer_parameters> const &params)
    {
      if (!has_floating_directions()) return;

      // construct singular vectors
      whole_deltas.clear();
      af::shared<scalar_t>
      scatterer_weights = weights(normal_eqns,
                  jacobian_transpose_matching_grad_fc,
                  params);
      SMTBX_ASSERT(params.size() == scatterer_weights.size())
                  (params.size())
                  (scatterer_weights.size());
      for (int k=0; k<deltas.size(); ++k) {
        af::shared<scalar_t> whole_delta(af::reserve(5*params.size()));
        std::back_insert_iterator< af::shared<scalar_t> > g(whole_delta);
        for (std::size_t i=0; i<params.size(); ++i) {
          scatterer_type const *sc = params[i].scatterer;
          constraints::index_range
          site_indices = params[i].site->component_indices_for(sc);
          BOOST_FOREACH (constraints::asu_parameter const *p,
                         params[i].ordered())
          {
            vec3_t v = scatterer_weights[i]*deltas[k];
            constraints::index_range
            indices = p->component_indices_for(sc);
            if (indices == site_indices) {
              std::copy(v.begin(), v.end(), g);
            }
            else {
              std::fill_n(g, indices.size(), scalar_t(0));
            }
          }
        }
        whole_deltas.push_back(whole_delta);
      }
      for (int i=0; i<whole_deltas.size(); ++i) {
        af::shared<scalar_t> reparametrised_delta =
        jacobian_transpose_matching_grad_fc*whole_deltas[i].const_ref();
        normal_eqns.add_equation(0, reparametrised_delta.const_ref(), 1);
      }
    }
  };

}}}


#endif // Guard
