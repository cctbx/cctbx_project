// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MATH_MIN_MAX_H
#define CCTBX_MATH_MIN_MAX_H

namespace cctbx { namespace math {

  template <typename NumType1, typename NumType2>
  NumType1&
  update_min(NumType1& m, NumType2 const& x)
  {
    if (m > x) m = x;
    return m;
  }

  template <typename NumType1, typename NumType2>
  NumType1&
  update_max(NumType1& m, NumType2 const& x)
  {
    if (m < x) m = x;
    return m;
  }

}} // namespace cctbx::math

#endif // CCTBX_MATH_MIN_MAX_H
