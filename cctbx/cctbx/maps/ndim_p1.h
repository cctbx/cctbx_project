// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_P1_ACCESS_H
#define CCTBX_MAPS_P1_ACCESS_H

#include <cctbx/ndim.h>

namespace cctbx {

//! Algorithms for the handling of 3-dimensional %maps.
namespace maps {

  template <std::size_t D, typename Index1dType = c_index_1d<D> >
  class dimension_p1 : public carray<int, D>
  {
    public:
      typedef carray<int, D> index_tuple_type;

      dimension_p1() {};
      dimension_p1(const carray<int, D>& N) {
        std::copy(N.begin(), N.end(), this->begin());
      }
      dimension_p1(std::size_t n0) {
        this->elems[0] = n0;
      }
      dimension_p1(std::size_t n0, std::size_t n1) {
        this->elems[0] = n0;
        this->elems[1] = n1;
      }
      dimension_p1(std::size_t n0, std::size_t n1, std::size_t n2) {
        this->elems[0] = n0;
        this->elems[1] = n1;
        this->elems[2] = n2;
      }

      std::size_t size1d() const { return cctbx::vector::product(*this); }

      template <typename IndexTupleType>
      carray<int, D>
      p1_I(const IndexTupleType& I) const {
        carray<int, D> result;
        for(std::size_t i=0;i<size();i++) {
          result[i] = I[i] % this->elems[i];
          if (result[i] < 0) result[i] += this->elems[i];
        }
        return result;
      }

      template <typename IndexTupleType>
      std::size_t operator()(const IndexTupleType& I) const {
        return Index1dType()(*this, p1_I(I));
      }

  };

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_P1_ACCESS_H
