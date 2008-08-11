#ifndef CCTBX_SGTBX_SMITH_NORMAL_FORM_H
#define CCTBX_SGTBX_SMITH_NORMAL_FORM_H

#include <scitbx/matrix/row_echelon.h>

namespace cctbx { namespace sgtbx {

  template <typename IntType>
  void
  smith_normal_form(
    scitbx::mat_ref<IntType>& m,
    scitbx::mat_ref<IntType> const& p,
    scitbx::mat_ref<IntType> const& q)
  {
    if (p.begin()) p.set_identity();
    if (q.begin()) q.set_identity();
    for (;;)
    {
      scitbx::matrix::row_echelon::form_t(m, p);
      if (m.is_diagonal()) break;
      m.transpose_in_place();

      scitbx::matrix::row_echelon::form_t(m, q);
      if (m.is_diagonal()) break;
      m.transpose_in_place();
    }
    if (q.begin()) q.transpose_square_in_place();
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SMITH_NORMAL_FORM_H
