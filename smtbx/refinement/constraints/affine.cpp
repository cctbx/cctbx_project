#include <smtbx/refinement/constraints/affine.h>

namespace smtbx { namespace refinement { namespace constraints {

  void affine_scalar_parameter::linearise(uctbx::unit_cell const &unit_cell,
                                          sparse_matrix_type *jacobian_transpose)
  {
    value = b;
    for(std::size_t i=0; i<n_arguments(); i++) value += a[i]*u(i)->value;
    if(!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for(std::size_t i=0; i<n_arguments(); i++) {
      jt.col(index()) += a[i]*jt.col(argument(i)->index());
    }
  }

}}}
