#include <smtbx/refinement/constraints/shared.h>

#include <boost/lambda/lambda.hpp>
#include <scitbx/array_family/tiny_algebra.h>

namespace smtbx { namespace refinement { namespace constraints {
  // u_star
  void
  shared_u_star
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    u_star_parameter *u = original();
    value = u->value;

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int j=0; j<6; ++j) {
      jt.col(index() + j) = jt.col(u->index() + j);
    }
  }

  // u_iso
  void
  shared_u_iso
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    scalar_parameter *u_iso = original();
    value = u_iso->value;

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.col(index()) = jt.col(u_iso->index());
  }

  // site
  void
  shared_site
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    site_parameter *site = original();
    value = site->value;

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int i=0; i < 3; i++)
      jt.col(index()+i) = jt.col(site->index()+i);
  }

}}}
