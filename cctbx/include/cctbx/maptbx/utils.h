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

  //! Miller index element corresponding to 1-dimensional array index.
  template <typename IntegerType>
  inline
  IntegerType
  ih_as_h(IntegerType ih, std::size_t n_real)
  {
    if (ih <= n_real/2) return ih;
    return ih - n_real;
  }

  //! 1-dimensional array index corresponding to Miller index element.
  /*! Returns -1 if h is out of range (see code).
   */
  template <typename IntegerType>
  inline
  IntegerType
  h_as_ih_exact(IntegerType h, IntegerType n_complex, bool positive_only)
  {
    if (positive_only) {
      if (0 > h || h >= n_complex) return -1;
    }
    else {
      IntegerType m = (n_complex - 1) / 2;
      if (-m > h || h > m) return -1;
      else if (h < 0) return h + n_complex;
    }
    return h;
  }

  //! 3-dimensional array indices corresponding to Miller index.
  /*! Result is -1 for out-of-range elements (see code).
   */
  template <typename IndexTypeN>
  af::int3
  h_as_ih_exact_array(bool anomalous_flag,
                      miller::index<> const& h,
                      IndexTypeN const& n_complex)
  {
    af::int3 ih;
    bool positive_only[] = {false, false, !anomalous_flag};
    for(std::size_t i=0;i<3;i++) {
      ih[i] = h_as_ih_exact(h[i], n_complex[i], positive_only[i]);
    }
    return ih;
  }

  //! 1-dimensional array index corresponding to Miller index element.
  /*! Applies modulus operation (see code).
   */
  template <typename IntegerType>
  inline
  IntegerType
  h_as_ih_mod(IntegerType h, IntegerType const& n_real)
  {
    h %= n_real;
    if (h < 0) return h + n_real;
    return h;
  }

  //! 3-dimensional array indices corresponding to Miller index.
  /*! Applies modulus operation (see code).
      <p>
      See also: structure_factors::to_map
   */
  template <typename IndexTypeN>
  IndexTypeN
  h_as_ih_mod_array(miller::index<> const& h, IndexTypeN const& n_real)
  {
    IndexTypeN ih;
    for(std::size_t i=0;i<3;i++) {
      ih[i] = h_as_ih_mod(h[i], n_real[i]);
    }
    return ih;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_UTILS_H
