#include <smtbx/refinement/constraints/special_position.h>

namespace smtbx { namespace refinement { namespace constraints {

  /**** Sites ****/

  void special_position_site::linearise(uctbx::unit_cell const &unit_cell,
                                        sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<3> const &p = independent_params();

    value = site_constraints.all_params(p.value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    af::const_ref<double, af::mat_grid>
    compact_jt = site_constraints.gradient_sum_matrix();
    for (int j=0; j<3; ++j) {
      for (int i=0; i<site_constraints.n_independent_params(); ++i) {
        jt(p.index() + i, index() + j) = compact_jt(i, j);
      }
    }
  }

  /**** ADP ****/

  void special_position_cartesian_adp
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<6> const &p = independent_params();

    value = adp_constraints.all_params(p.value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jac_t = *jacobian_transpose;
    sparse_matrix_type const &jac = adp_constraints.jacobian();
    typedef sparse_matrix_type::const_row_iterator iter_t;
    for (int j=0; j<jac.n_cols(); ++j) {
      for (iter_t i=jac.col(j).begin(); i != jac.col(j).end(); ++i) {
        jac_t(j, index() + i.index()) = *i;
      }
    }
  }

}}}
