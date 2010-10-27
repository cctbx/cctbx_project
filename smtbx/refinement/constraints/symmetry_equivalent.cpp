#include <smtbx/refinement/constraints/symmetry_equivalent.h>

namespace smtbx { namespace refinement { namespace constraints {

  // Symmetry equivalent site

  void symmetry_equivalent_site_parameter
  ::linearise(uctbx::unit_cell const &unit_cell,
              sparse_matrix_type *jacobian_transpose)
  {
    site_parameter *site = original();

    value = op*site->value;

    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.assign_block(local_jt, site->index(), index());
  }

}}}
