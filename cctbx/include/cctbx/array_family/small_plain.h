// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
#define CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H

#include <cctbx/array_family/tiny.h>
#include <cctbx/array_family/auto_allocator.h>

namespace cctbx { namespace af {

  // Automatic allocation, variable size.
  template <typename ElementType, std::size_t N>
  class small_plain
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      small_plain()
        : m_size(0)
      {}

      explicit
      small_plain(size_type const& sz)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
        std::uninitialized_fill_n(begin(), sz, ElementType());
        m_size = sz;
      }

      // non-std
      small_plain(size_type const& sz, reserve_flag)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
      }

      small_plain(size_type const& sz, ElementType const& x)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
        std::uninitialized_fill_n(begin(), sz, x);
        m_size = sz;
      }

      // non-std
      template <typename FunctorType>
      small_plain(size_type const& sz, init_functor<FunctorType> const& ftor)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
        ftor.held(begin(), sz);
        m_size = sz;
      }

      small_plain(const ElementType* first, const ElementType* last)
        : m_size(0)
      {
        if (N < last - first) throw_range_error();
        std::uninitialized_copy(first, last, begin());
        m_size = last - first;
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      small_plain(const OtherElementType* first, const OtherElementType* last)
        : m_size(0)
      {
        if (N < last - first) throw_range_error();
        uninitialized_copy_typeconv(first, last, begin());
        m_size = last - first;
      }
#endif

      small_plain(small_plain<ElementType, N> const& other)
        : m_size(0)
      {
        std::uninitialized_copy(other.begin(), other.end(), begin());
        m_size = other.m_size;
      }

      template <typename OtherArrayType>
      small_plain(array_adaptor<OtherArrayType> const& a_a)
        : m_size(0)
      {
        OtherArrayType const& a = *(a_a.pointee);
        if (a.size() > N) throw_range_error();
        for(std::size_t i=0;i<a.size();i++) push_back(a[i]);
      }

      ~small_plain() { clear(); }

      small_plain<ElementType, N>&
      operator=(small_plain<ElementType, N> const& other)
      {
        clear();
        std::uninitialized_copy(other.begin(), other.end(), begin());
        m_size = other.m_size;
        return *this;
      }

      size_type size() const { return m_size; }
      bool empty() const { if (size() == 0) return true; return false; }
      static size_type max_size() { return N; }
      static size_type capacity() { return N; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(small_plain,
        ((ElementType*)(m_elems.buffer)), m_size) // fix this

      CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), m_size)

      void swap(small_plain<ElementType, N>& other) {
        std::swap(*this, other);
      }

      void reserve(size_type const& sz) {
        if (N < sz) throw_range_error();
      }

#     include <cctbx/array_family/push_back_etc.h>

    protected:
      void m_insert_overflow(ElementType* pos,
                             size_type const& n, ElementType const& x,
                             bool at_end) {
        throw_range_error();
      }

      void m_insert_overflow(ElementType* pos,
                             const ElementType* first,
                             const ElementType* last) {
        throw_range_error();
      }

      void m_set_size(size_type const& sz) {
        m_size = sz;
      }

      void m_incr_size(size_type const& n) {
        m_size += n;
      }

      void m_decr_size(size_type const& n) {
        m_size -= n;
      }

      detail::auto_allocator<ElementType, N> m_elems;
      size_type m_size;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
