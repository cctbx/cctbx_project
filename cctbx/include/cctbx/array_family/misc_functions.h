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
  approx_equal_scaled(FloatType const& a,
                      FloatType const& b,
                      FloatType const& scaled_tolerance) {
    FloatType diff = a - b;
    if (diff < 0.) diff = -diff;
    if (diff < scaled_tolerance) return true;
    return false;
  }

  //! Test if 2*abs((a-b)/(a+b)) < tolerance.
  template <class FloatType>
  bool
  approx_equal_unscaled(FloatType const& a,
                        FloatType const& b,
                        FloatType const& tolerance) {
    FloatType sum = a + b;
    cctbx_assert(sum != 0);
    FloatType diff = a - b;
    FloatType ratio = diff / sum;
    if (ratio < 0) ratio = -ratio;
    if (FloatType(2) * ratio < tolerance) return true;
    return false;
  }

  //! Helper function object for array operations.
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2,
            typename ArgumentType3>
  struct functor_approx_equal_scaled {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType1 const& x,
                          ArgumentType2 const& y,
                          ArgumentType3 const& z) const {
    return ResultType(approx_equal_scaled(x, y, z)); }
  };

  //! Helper function object for array operations.
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2,
            typename ArgumentType3>
  struct functor_approx_equal_unscaled {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType1 const& x,
                          ArgumentType2 const& y,
                          ArgumentType3 const& z) const {
    return ResultType(approx_equal_unscaled(x, y, z)); }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_MISC_FUNCTIONS_H
