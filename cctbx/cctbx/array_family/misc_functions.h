// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: moved from utils.h to array_family/misc_functions.h (rwgk)
     2001 Oct 12: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_MISC_FUNCTIONS_H
#define CCTBX_ARRAY_FAMILY_MISC_FUNCTIONS_H

namespace cctbx { namespace af {

  //! Test if abs(a-b) < scaled_tolerance.
  template <class FloatType>
  bool
  approx_equal_scaled(const FloatType& a,
                      const FloatType& b,
                      const FloatType& scaled_tolerance) {
    FloatType diff = a - b;
    if (diff < 0.) diff = -diff;
    if (diff < scaled_tolerance) return true;
    return false;
  }

  //! Test if 2*abs((a-b)/(a+b)) < tolerance.
  template <class FloatType>
  bool
  approx_equal_unscaled(const FloatType& a,
                        const FloatType& b,
                        const FloatType& tolerance) {
    FloatType sum = a + b;
    cctbx_assert(sum != 0);
    FloatType diff = a - b;
    FloatType ratio = diff / sum;
    if (ratio < 0) ratio = -ratio;
    if (FloatType(2) * ratio < tolerance) return true;
    return false;
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_MISC_FUNCTIONS_H
