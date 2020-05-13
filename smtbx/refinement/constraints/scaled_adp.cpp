#include <smtbx/refinement/constraints/scaled_adp.h>

namespace smtbx { namespace refinement { namespace constraints {

void
scalar_scaled_u_star_parameter
::linearise(uctbx::unit_cell const &unit_cell,
            sparse_matrix_type *jacobian_transpose)
{
  const scitbx::sym_mat3<double> *ref = reference();
  independent_scalar_parameter *scalar_ = scalar();
  value = *ref * scalar_->value;

  if (!jacobian_transpose) return;
  sparse_matrix_type &jt = *jacobian_transpose;
  std::size_t i = scalar_->index();
  for (int k=0; k<6; ++k) {
    jt(i, index() + k) = (*ref)[k];
  }
}

void
scalar_scaled_u_iso_parameter
::linearise(uctbx::unit_cell const &unit_cell,
            sparse_matrix_type *jacobian_transpose)
{
  const double ref = reference();
  independent_scalar_parameter *scalar_ = scalar();
  value = ref * scalar_->value;

  if (!jacobian_transpose) return;
  sparse_matrix_type &jt = *jacobian_transpose;
  std::size_t i = scalar_->index();
  jt(i, index()) = ref;
}
}}}
