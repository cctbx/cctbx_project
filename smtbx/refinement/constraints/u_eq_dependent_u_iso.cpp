#include <smtbx/refinement/constraints/u_eq_dependent_u_iso.h>

#include <boost/lambda/lambda.hpp>
#include <scitbx/array_family/tiny_algebra.h>

namespace smtbx { namespace refinement { namespace constraints {

  void
  u_iso_proportional_to_pivot_u_eq
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    af::double6 f = unit_cell.u_star_to_u_iso_linear_form();
    f *= multiplier;
    u_star_parameter *u = pivot_u();
    value = f * u->value;

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int j=0; j<6; ++j) {
      jt.col(index()) = jt.col(index()) + f[j]*jt.col(u->index() + j);
    }
  }

}}}
