#include <smtbx/refinement/constraints/occupancy.h>

#include <boost/lambda/lambda.hpp>
#include <scitbx/array_family/tiny_algebra.h>

namespace smtbx { namespace refinement { namespace constraints {

  void
  dependent_occupancy
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    scalar_parameter *o = reference();
    if( as_one )
      value = multiplier*o->value;
    else
      value = multiplier*(1.0 - o->value);

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    if( as_one )
      jt.col(index()) = multiplier*jt.col(o->index());
    else
      jt.col(index()) = (-multiplier)*jt.col(o->index());
  }

}}}
