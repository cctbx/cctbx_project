// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 12: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_UTILS_H
#define CCTBX_UTILS_H

namespace cctbx {

  //! Test if abs(a-b) < scaled_tolerance.
  template <class FloatType>
  bool
  approx_equal(const FloatType& a,
               const FloatType& b,
               const FloatType& scaled_tolerance) {
    FloatType diff = a - b;
    if (diff < 0.) diff = -diff;
    if (diff < scaled_tolerance) return true;
    return false;
  }

}

#endif // CCTBX_UTILS_H
