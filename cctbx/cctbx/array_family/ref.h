// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_REF_H
#define CCTBX_ARRAY_FAMILY_REF_H

#include <stdexcept>

// FIXES for broken compilers
#include <boost/config.hpp>

// XXX
#include <iostream>
#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

#include <cctbx/array_family/grid_accessor.h>
#include <cctbx/array_family/ref_helpers.h>
#include <cctbx/array_family/versa_helpers.h>

namespace cctbx { namespace af {

  inline
  void throw_range_error() {
    throw std::range_error("array_family");
  }

  template <typename InputElementType,
            typename OutputElementType>
  OutputElementType*
  copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
    OutputElementType* p = result;
    while (first != last) *p++ = OutputElementType(*first++);
    return result;
  }

  template <typename ElementType,
            typename AccessorType = grid_accessor<1> >
  class const_ref
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;

      const_ref()
        : m_begin(0) {
        m_accessor.init_default();
      }
      const_ref(const ElementType* begin, accessor_type ac)
        : m_begin(begin), m_accessor(ac)
      {}
      const_ref(const void* begin, accessor_type ac)
        : m_begin(reinterpret_cast<const ElementType*>(begin)),
          m_accessor(ac)
      {}
      // convenience constructors
      const_ref(
        const ElementType* begin, size_type n0)
        : m_begin(begin), m_accessor(n0)
      {}
      const_ref(
        const ElementType* begin, size_type n0, size_type n1)
        : m_begin(begin), m_accessor(n0, n1)
      {}
      const_ref(
        const ElementType* begin, size_type n0, size_type n1, size_type n2)
        : m_begin(begin), m_accessor(n0, n1, n2)
      {}

      const accessor_type& accessor() { return m_accessor; }
      size_type size() const { return m_accessor.size1d(); }

      const ElementType* begin() const { return m_begin; }
      const ElementType* end() const { return m_begin + size(); }
      const ElementType& front() const { return m_begin[0]; }
      const ElementType& back() const { return m_begin[size()-1]; }

      const ElementType& operator[](size_type i) const { return m_begin[i]; }

      const ElementType& at(size_type i) const {
        if (i >= size()) throw_range_error();
        return m_begin[i];
      }

      const void* handle() const {
        return reinterpret_cast<const void*>(m_begin);
      }

      const value_type& operator()(const index_type& i) const {
        return this->begin()[m_accessor(i)];
      }

      // Convenience operator()

      const value_type& operator()(int i0) const {
        return operator()(index_type(i0));
      }
      const value_type& operator()(int i0,
                                   int i1) const {
        return operator()(index_type(i0, i1));
      }
      const value_type& operator()(int i0,
                                   int i1,
                                   int i2) const {
        return operator()(index_type(i0, i1, i2));
      }

    private:
      const ElementType* m_begin;
      accessor_type m_accessor;
  };

  template <typename ElementType,
            typename AccessorType = grid_accessor<1> >
  class ref
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;

      ref()
        : m_begin(0) {
        m_accessor.init_default();
      }
      ref(ElementType* begin, accessor_type ac)
        : m_begin(begin), m_accessor(ac)
      {}
      ref(void* begin, accessor_type ac)
        : m_begin(reinterpret_cast<ElementType*>(begin)),
          m_accessor(ac)
      {}
      // convenience constructors
      ref(
        ElementType* begin, size_type n0)
        : m_begin(begin), m_accessor(n0)
      {}
      ref(
        ElementType* begin, size_type n0, size_type n1)
        : m_begin(begin), m_accessor(n0, n1)
      {}
      ref(
        ElementType* begin, size_type n0, size_type n1, size_type n2)
        : m_begin(begin), m_accessor(n0, n1, n2)
      {}

      const accessor_type& accessor() { return m_accessor; }
      size_type size() const { return m_accessor.size1d(); }

      ElementType* begin() const { return m_begin; }
      ElementType* end() const { return m_begin + size(); }
      ElementType& front() const { return m_begin[0]; }
      ElementType& back() const { return m_begin[size()-1]; }

      ElementType& operator[](size_type i) const { return m_begin[i]; }

      ElementType& at(size_type i) const {
        if (i >= size()) throw_range_error();
        return m_begin[i];
      }

      void* handle() const {
        return reinterpret_cast<void*>(m_begin);
      }

      operator const_ref<ElementType, AccessorType>() const {
        return const_ref<ElementType, AccessorType>(m_begin, m_accessor);
      }

            value_type& operator()(const index_type& i) const {
        return this->begin()[m_accessor(i)];
      }

      // Convenience operator()

            value_type& operator()(int i0) const {
        return operator()(index_type(i0));
      }
            value_type& operator()(int i0,
                                   int i1) const {
        return operator()(index_type(i0, i1));
      }
            value_type& operator()(int i0,
                                   int i1,
                                   int i2) const {
        return operator()(index_type(i0, i1, i2));
      }

    private:
      ElementType* m_begin;
      accessor_type m_accessor;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_REF_H
