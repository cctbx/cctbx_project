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

#include <algorithm>
#include <cctbx/vecref.h>
#include <cctbx/vector/reductions.h>

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

  template <std::size_t D, typename Index1dType = c_index_1d<D> >
  class dimension : public carray<int, D>
  {
    public:
      dimension() {};
      dimension(const carray<int, D>& N) {
        std::copy(N.begin(), N.end(), this->begin());
      }
      dimension(std::size_t n0) {
        this->elems[0] = n0;
      }
      dimension(std::size_t n0, std::size_t n1) {
        this->elems[0] = n0;
        this->elems[1] = n1;
      }
      dimension(std::size_t n0, std::size_t n1, std::size_t n2) {
        this->elems[0] = n0;
        this->elems[1] = n1;
        this->elems[2] = n2;
      }

      std::size_t size1d() const { return cctbx::vector::product(*this); }

      template <typename IndexTuple>
      std::size_t operator()(const IndexTuple& I) const {
        return Index1dType()(*this, I);
      }

      template <typename IndexTuple>
      bool is_valid_index(const IndexTuple& I) const {
        if (I.size() != size()) return false;
        for(std::size_t j=0;j<size();j++) {
          std::size_t i = I[j];
          if (i >= this->elems[j]) return false;
        }
        return true;
      }
  };

  template <typename ValueType, typename DimensionType>
  class vecrefnd : public vecref<ValueType>
  {
    public:
      typedef ValueType value_type;
      typedef DimensionType dimension_type;

      vecrefnd() {}
      template <typename IteratorOrPointerType>
      vecrefnd(IteratorOrPointerType begin, const dimension_type& dim)
        : vecref<ValueType>(begin, dim.size1d()), m_dim(dim)
      {}
      vecrefnd(void* begin, const dimension_type& dim)
        : vecref<ValueType>(begin, dim.size1d()), m_dim(dim)
      {}

      const dimension_type& dim() const { return m_dim; }

      template <typename IndexTuple>
      value_type& operator()(const IndexTuple& I) const {
        return m_begin[m_dim(I)];
      }

    protected:
      dimension_type m_dim;
  };

}

#endif // CCTBX_NDIM_H
