/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Generalized cctbx/sgtbx/mod.h (rwgk)
     2002 Sep: Fragment from cctbx/sgtbx/math.h (rwgk)
 */

#ifndef CCTBX_MATH_MOD_H
#define CCTBX_MATH_MOD_H

namespace cctbx { namespace math {

  template <typename IntType>
  inline
  IntType
  mod_positive(IntType ix, IntType const& iy)
  {
    if (iy > 0) {
      ix %= iy;
      if (ix < 0) ix += iy;
    }
    return ix;
  }

  template <typename IntType>
  inline
  IntType
  mod_short(IntType ix, IntType const& iy)
  {
        ix = mod_positive(ix, iy);
    if (ix > iy / 2)
        ix -= iy;
    return ix;
  }

}} // namespace cctbx::math

#endif // CCTBX_MATH_MOD_H
