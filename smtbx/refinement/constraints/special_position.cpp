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
    local_jt = site_constraints.gradient_sum_matrix();
    for (int j=0; j<3; ++j) {
      for (int i=0; i<site_constraints.n_independent_params(); ++i) {
        if (local_jt(i, j)) jt(p.index() + i, index() + j) = local_jt(i, j);
      }
    }
  }

  /**** ADP ****/

  void special_position_u_star_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<6> const &p = independent_params();

    value = adp_constraints.all_params(p.value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    af::const_ref<double, af::mat_grid>
    local_jt = adp_constraints.gradient_sum_matrix();
    for (int j=0; j<6; ++j) {
      for (int i=0; i<adp_constraints.n_independent_params(); ++i) {
        if (local_jt(i, j)) jt(p.index() + i, index() + j) = local_jt(i, j);
      }
    }
  }

}}}
