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

#include <algorithm>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/auto_allocator.h>
#include <cctbx/array_family/type_traits.h>

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
      small_plain(const size_type& sz)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
        std::uninitialized_fill_n(begin(), sz, ElementType());
        m_size = sz;
      }

      // non-std
      small_plain(const size_type& sz, reserve_flag)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
      }

      // non-std
      small_plain(const size_type& sz, no_initialization_flag)
        : m_size(0)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        if (N < sz) throw_range_error();
        m_size = sz;
      }

      small_plain(const size_type& sz, const ElementType& x)
        : m_size(0)
      {
        if (N < sz) throw_range_error();
        std::uninitialized_fill_n(begin(), sz, x);
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

      small_plain(const small_plain<ElementType, N>& other)
        : m_size(0)
      {
        std::uninitialized_copy(other.begin(), other.end(), begin());
        m_size = other.m_size;
      }

      ~small_plain() { clear(); }

      small_plain<ElementType, N>&
      operator=(const small_plain<ElementType, N>& other)
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

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
        ((ElementType*)(m_elems.buffer)), m_size) // fix this

      CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), N)

      void swap(small_plain<ElementType, N>& other) {
        std::swap(*this, other);
      }

      void reserve(const size_type& sz) {
        if (N < sz) throw_range_error();
      }

      // non-std
      void set_size_back_door(const size_type& sz) {
        m_set_size(sz);
      }

#     include <cctbx/array_family/push_back_etc.h>

    protected:
      void m_insert_overflow(ElementType* pos,
                             const size_type& n, const ElementType& x,
                             bool at_end) {
        throw_range_error();
      }

      void m_insert_overflow(ElementType* pos,
                             const ElementType* first,
                             const ElementType* last) {
        throw_range_error();
      }

      void m_set_size(const size_type& sz) {
        m_size = sz;
      }

      void m_incr_size(const size_type& n) {
        m_size += n;
      }

      void m_decr_size(const size_type& n) {
        m_size -= n;
      }

      detail::auto_allocator<ElementType, N> m_elems;
      size_type m_size;
  };

}} // namespace cctbx::af

#include <cctbx/array_family/small_plain_apply.h>

#endif // CCTBX_ARRAY_FAMILY_SMALL_PLAIN_H
