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
#include <cctbx/carray.h>
#include <cctbx/shared_storage.h>
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
      typedef carray<int, D> index_tuple_type;

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

      template <typename IndexTupleType>
      std::size_t operator()(const IndexTupleType& I) const {
        return Index1dType()(*this, I);
      }

      template <typename IndexTupleType>
      bool is_valid_index(const IndexTupleType& I) const {
        if (I.size() != this->size()) return false;
        for(std::size_t j=0;j<this->size();j++) {
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
      vecrefnd(char* begin, const dimension_type& dim)
        : vecref<ValueType>(begin, dim.size1d()), m_dim(dim)
      {}

      const dimension_type& dim() const { return m_dim; }

      const vecref<ValueType>& as_1d() const { return *this; }

      template <typename IndexObjectType>
      value_type& operator()(const IndexObjectType& I) const {
        return this->m_begin[this->m_dim(I)];
      }

      // Convenience constructors

      template <typename IteratorOrPointerType>
      vecrefnd(IteratorOrPointerType begin,
               std::size_t n0)
        : vecref<ValueType>(begin, n0), m_dim(n0)
      {}
      template <typename IteratorOrPointerType>
      vecrefnd(IteratorOrPointerType begin,
               std::size_t n0,
               std::size_t n1)
        : vecref<ValueType>(begin, n0*n1), m_dim(n0, n1)
      {}
      template <typename IteratorOrPointerType>
      vecrefnd(IteratorOrPointerType begin,
               std::size_t n0,
               std::size_t n1,
               std::size_t n2)
        : vecref<ValueType>(begin, n0*n1*n2), m_dim(n0, n1, n2)
      {}

      // Convenience operator()

      value_type& operator()(std::size_t i0) const {
        return this->m_begin[this->m_dim(carray<int, 1>(i0))];
      }
      value_type& operator()(std::size_t i0,
                             std::size_t i1) const {
        return this->m_begin[this->m_dim(carray<int, 2>(i0, i1))];
      }
      value_type& operator()(std::size_t i0,
                             std::size_t i1,
                             std::size_t i2) const {
        return this->m_begin[this->m_dim(carray<int, 3>(i0, i1, i2))];
      }

    protected:
      dimension_type m_dim;
  };

  template <typename ValueType, typename DimensionType>
  class shared_storage_nd : public shared_storage<ValueType>
  {
    public:
      typedef ValueType value_type;
      typedef DimensionType dimension_type;
      typedef typename shared_storage<value_type>::handle_type handle_type;

      shared_storage_nd() {}
      shared_storage_nd(const dimension_type& dim)
        : shared_storage<value_type>(dim.size1d()), m_dim(dim)
      {}
      shared_storage_nd(const handle_type& handle, const dimension_type& dim)
        : shared_storage<value_type>(handle, dim.size1d()), m_dim(dim)
      {}

      const dimension_type& dim() const { return m_dim; }

      vecrefnd<      value_type, dimension_type> ref()       {
        return vecrefnd<      value_type, dimension_type>(
          this->m_begin, this->m_dim);
      }
      vecrefnd<const value_type, dimension_type> ref() const {
        return vecrefnd<const value_type, dimension_type>(
          this->m_begin, this->m_dim);
      }

            shared_storage<value_type>& as_1d()       { return *this; }
      const shared_storage<value_type>& as_1d() const { return *this; }

      template <typename IndexObjectType>
      value_type& operator()(const IndexObjectType& I) const {
        return this->m_begin[this->m_dim(I)];
      }

      // Convenience constructors

      shared_storage_nd(std::size_t n0)
        : shared_storage<value_type>(n0), m_dim(n0)
      {}
      shared_storage_nd(std::size_t n0,
                        std::size_t n1)
        : shared_storage<value_type>(n0*n1), m_dim(n0, n1)
      {}
      shared_storage_nd(std::size_t n0,
                        std::size_t n1,
                        std::size_t n2)
        : shared_storage<value_type>(n0*n1*n2), m_dim(n0, n1, n2)
      {}

      // Convenience operator()

      value_type& operator()(std::size_t i0) const {
        return this->m_begin[this->m_dim(carray<int, 1>(i0))];
      }
      value_type& operator()(std::size_t i0,
                             std::size_t i1) const {
        return this->m_begin[this->m_dim(carray<int, 2>(i0, i1))];
      }
      value_type& operator()(std::size_t i0,
                             std::size_t i1,
                             std::size_t i2) const {
        return this->m_begin[this->m_dim(carray<int, 3>(i0, i1, i2))];
      }

    protected:
      dimension_type m_dim;
  };

}

#endif // CCTBX_NDIM_H
