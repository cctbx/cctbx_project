// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_NDIM_H
#define CCTBX_NDIM_H

namespace cctbx {

  template <std::size_t N>
  struct c_index_1d {
    template <typename ExtendArrayType, typename IndexArrayType>
    std::size_t operator()(const ExtendArrayType& e, const IndexArrayType& i) {
      return c_index_1d<N-1>()(e, i) * e[N-1] + i[N-1];
    }
  };

  template<>
  struct c_index_1d<1> {
    template <typename ExtendArrayType, typename IndexArrayType>
    std::size_t operator()(const ExtendArrayType& e, const IndexArrayType& i) {
      return i[0];
    }
  };

  template <std::size_t N>
  struct fortran_index_1d {
    template <typename ExtendArrayType, typename IndexArrayType>
    std::size_t operator()(const ExtendArrayType& e, const IndexArrayType& i) {
      return fortran_index_1d<N-1>()(e, i) * e[e.size()-N] + i[e.size()-N];
    }
  };

  template<>
  struct fortran_index_1d<1> {
    template <typename ExtendArrayType, typename IndexArrayType>
    std::size_t operator()(const ExtendArrayType& e, const IndexArrayType& i) {
      return i[e.size()-1];
    }
  };

}

#endif // CCTBX_NDIM_H
