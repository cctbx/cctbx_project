/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored part of sgtbx/math.cpp (rwgk)
 */

#ifndef CCTBX_SGTBX_SMITH_NORMAL_FORM_H
#define CCTBX_SGTBX_SMITH_NORMAL_FORM_H

#include <cctbx/sgtbx/row_echelon.h>

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
      row_echelon::form_t(m, p);
      if (m.is_diagonal()) break;
      m.transpose_in_place();

      row_echelon::form_t(m, q);
      if (m.is_diagonal()) break;
      m.transpose_in_place();
    }
    if (q.begin()) q.transpose_square_in_place();
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SMITH_NORMAL_FORM_H
