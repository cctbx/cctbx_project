/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Fragment from cctbx/maps/accessors.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_UTILS_H
#define CCTBX_MAPTBX_UTILS_H

#include <cstddef>

namespace cctbx { namespace maptbx {

  //! 1-dimensional array index corresponding to Miller index element.
  template <typename IntegerType>
  inline
  IntegerType
  h_as_ih(IntegerType h, std::size_t n, bool positive_only)
  {
    if (positive_only) {
      if (0 > h || h >= n) return -1;
    }
    else {
      IntegerType m = (n - 1) / 2;
      if (-m > h || h > m) return -1;
      else if (h < 0) return h + n;
    }
    return h;
  }

  template <typename IndexTypeN>
  miller::index<>
  h_as_ih_array(bool anomalous_flag,
                miller::index<> const& h,
                IndexTypeN const& n)
  {
    miller::index<> ih;
    bool positive_only[] = {false, false, !anomalous_flag};
    for(std::size_t i=0;i<3;i++) {
      ih[i] = h_as_ih(h[i], n[i], positive_only[i]);
    }
    return ih;
  }

  //! Miller index element corresponding to 1-dimensional array index.
  template <typename IntegerType>
  inline
  IntegerType
  ih_as_h(IntegerType ih, std::size_t n)
  {
    if (ih <= n/2) return ih;
    return ih - n;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_UTILS_H
