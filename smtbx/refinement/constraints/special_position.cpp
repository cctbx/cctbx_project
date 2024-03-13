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

  /**** Anharmonic ADP ****/

  void special_position_anharmonic_adp_parameter
    ::linearise(uctbx::unit_cell const &unit_cell,
      sparse_matrix_type *jacobian_transpose)
  {
    int order = this->scatterer->anharmonic_adp->order;
    independent_vector_parameter const &p = independent_params();
    af::shared<double> tmp(c_count);
    for (size_t i = 0; i < c_count; i++) {
      tmp[i] = p.value[i];
    }
    tmp = tensor_r3_constraints.all_params(tmp);
    for (size_t i = 0; i < 10; i++) {
      value[i] = tmp[i];
    }
    if (order > 3) {
      tmp.resize(p.size() - c_count);
      for (size_t i = c_count; i < p.size(); i++) {
        tmp[i - c_count] = p.value[i];
      }
      tmp = tensor_r4_constraints.all_params(tmp);
      for (size_t i = 0; i < 15; i++) {
        value[i + 10] = tmp[i];
      }
    }

    if (!jacobian_transpose) {
      return;
    }
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.assign_block(tensor_r3_constraints.gradient_sum_matrix(), p.index(), index());
    if (order > 3) {
      jt.assign_block(tensor_r4_constraints.gradient_sum_matrix(),
        p.index() + c_count, index() + 10);
    }
  }

}}}
