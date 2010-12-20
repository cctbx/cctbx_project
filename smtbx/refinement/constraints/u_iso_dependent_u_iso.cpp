#include <smtbx/refinement/constraints/u_iso_dependent_u_iso.h>

#include <boost/lambda/lambda.hpp>
#include <scitbx/array_family/tiny_algebra.h>

namespace smtbx { namespace refinement { namespace constraints {

  void
  u_iso_proportional_to_pivot_u_iso
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    scalar_parameter *u_iso = pivot_u_iso();
    value = multiplier * u_iso->value;

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.col(index()) = jt.col(index()) + multiplier*jt.col(u_iso->index());
  }

}}}
