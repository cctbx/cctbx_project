#include <smtbx/refinement/constraints/special_position.h>

namespace smtbx { namespace refinement { namespace constraints {

  /**** Sites ****/

  void special_position_site_parameter::linearise(uctbx::unit_cell const &unit_cell,
                                        sparse_matrix_type *jacobian_transpose)
  {
    independent_small_vector_parameter<3> const &p = independent_params();

    value = site_constraints.all_params(p.value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.assign_block(site_constraints.gradient_sum_matrix(), p.index(), index());
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
    jt.assign_block(adp_constraints.gradient_sum_matrix(), p.index(), index());
  }

}}}
