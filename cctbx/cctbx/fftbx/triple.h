// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 03: fftbx started, based on fftpack41 (rwgk)
 */

#ifndef CCTBX_FFTBX_TRIPLE_H
#define CCTBX_FFTBX_TRIPLE_H

#include <boost/array.hpp>
#include <algorithm>

namespace cctbx { namespace fftbx {

  //! Triple of integers for 3-dimensional array sizes and indices.
  class triple : public boost::array<std::size_t, 3>
  {
    public:
      //! Default constructor. Triple elements are not initialized.
      triple() {};
      //! Initialization with three integers.
      triple(std::size_t ix, std::size_t iy, std::size_t iz) {
        elems[0] = ix;
        elems[1] = iy;
        elems[2] = iz;
      }
      //! Initialization with any boost::array of size 3.
      template <class T>
      triple(const boost::array<T, 3>& a) {
        for(std::size_t i=0;i<3;i++) {
          elems[i] = a[i];
        }
      }
      //! Determination of the largest element in the %triple.
      std::size_t max() const {
        return *std::max_element(begin(), end());
      }
      //! Product of the elements in the %triple.
      std::size_t product() const {
        return elems[0] * elems[1] * elems[2];
      }
  };

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_TRIPLE_H
