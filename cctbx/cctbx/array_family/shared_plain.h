// $Id$
/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: adaption of stlport code (rwgk)
 */
/* This code is derived in part from:
 *   STLport-4.5.3/stlport/stl/_vector.h
 *   STLport-4.5.3/stlport/stl/_vector.c
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Copyright (c) 1997
 * Moscow Center for SPARC Technology
 *
 * Copyright (c) 1999
 * Boris Fomitchev
 *
 * This material is provided "as is", with absolutely no warranty expressed
 * or implied. Any use is at your own risk.
 *
 * Permission to use or copy this software for any purpose is hereby granted
 * without fee, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is granted,
 * provided the above notices are retained, and a notice that the code was
 * modified is included with the above copyright notice.
 *
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H
#define CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H

#include <algorithm>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/type_traits.h>

namespace cctbx { namespace af {

  struct weak_ref_flag {};

  namespace detail {
    const std::size_t global_max_size(std::size_t(-1));
  }

  template <typename ElementType>
  class shared_plain
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      static size_type element_size() { return sizeof(ElementType); }

    protected:
      class handle_type {
        public:
          handle_type()
            : use_count(1), weak_count(0), size(0), capacity(0),
              data(0)
          {}

          handle_type(weak_ref_flag)
            : use_count(0), weak_count(1), size(0), capacity(0),
              data(0)
          {}

          explicit
          handle_type(const size_type& sz)
            : use_count(1), weak_count(0), size(0), capacity(sz),
              data(new char[sz])
          {}

          ~handle_type() {
            delete[] data;
          }

          void deallocate() {
            delete[] data;
            capacity = 0;
            data = 0;
          }

          void swap(handle_type& other) {
            std::swap(size, other.size);
            std::swap(capacity, other.capacity);
            std::swap(data, other.data);
          }

          size_type use_count;
          size_type weak_count;
          size_type size;
          size_type capacity;
          char* data;

        private:
          handle_type(const handle_type&);
          handle_type& operator=(const handle_type&);
      };

    public:
      shared_plain()
        : m_is_weak_ref(false),
          m_handle(new handle_type)
      {}

      explicit
      shared_plain(const size_type& sz)
        : m_is_weak_ref(false),
          m_handle(new handle_type(sz * element_size()))
      {
        std::uninitialized_fill_n(begin(), sz, ElementType());
        m_handle->size = m_handle->capacity;
      }

      // non-std
      shared_plain(const size_type& sz, reserve_flag)
        : m_is_weak_ref(false),
          m_handle(new handle_type(sz * element_size()))
      {}

      // non-std
      shared_plain(const size_type& sz, no_initialization_flag)
        : m_is_weak_ref(false),
          m_handle(new handle_type(sz * element_size()))
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->size = m_handle->capacity;
      }

      shared_plain(const size_type& sz, const ElementType& x)
        : m_is_weak_ref(false),
          m_handle(new handle_type(sz * element_size()))
      {
        std::uninitialized_fill_n(begin(), sz, x);
        m_handle->size = m_handle->capacity;
      }

      shared_plain(const ElementType* first, const ElementType* last)
        : m_is_weak_ref(false),
          m_handle(new handle_type((last - first) * element_size()))
      {
        std::uninitialized_copy(first, last, begin());
        m_handle->size = m_handle->capacity;
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      shared_plain(const OtherElementType* first, const OtherElementType* last)
        : m_is_weak_ref(false),
          m_handle(new handle_type((last - first) * element_size()))
      {
        uninitialized_copy_typeconv(first, last, begin());
        m_handle->size = m_handle->capacity;
      }
#endif

      // non-std: shallow copy semantics
      shared_plain(const shared_plain<ElementType>& other)
        : m_is_weak_ref(other.m_is_weak_ref),
          m_handle(other.m_handle)
      {
        if (m_is_weak_ref) m_handle->weak_count++;
        else               m_handle->use_count++;
      }

      // non-std: shallow copy semantics, weak reference
      shared_plain(const shared_plain<ElementType>& other, weak_ref_flag)
        : m_is_weak_ref(true),
          m_handle(other.m_handle)
      {
        m_handle->weak_count++;
      }

      // non-std
      explicit
      shared_plain(handle_type* other_handle)
        : m_is_weak_ref(false),
          m_handle(other_handle)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->use_count++;
      }

      // non-std
      explicit
      shared_plain(handle_type* other_handle, weak_ref_flag)
        : m_is_weak_ref(true),
          m_handle(other_handle)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->weak_count++;
      }

      ~shared_plain() {
        m_dispose();
      }

      // non-std: shallow copy semantics
      shared_plain<ElementType>&
      operator=(const shared_plain<ElementType>& other)
      {
        if (m_handle != other.m_handle) {
          m_dispose();
          m_is_weak_ref = other.m_is_weak_ref;
          m_handle = other.m_handle;
          if (m_is_weak_ref) m_handle->weak_count++;
          else               m_handle->use_count++;
        }
        return *this;
      }

      // non-std, const correctness is lost
      handle_type* handle() const {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        return m_handle;
      }

      size_type size() const { return m_handle->size / element_size(); }

      bool empty() const { if (size() == 0) return true; return false; }

      size_type capacity() const {
        return m_handle->capacity / element_size();
      }

      static size_type max_size() {
        return detail::global_max_size / element_size();
      }

      // non-std
      size_type use_count() const { return m_handle->use_count; }

      // non-std
      size_type weak_count() const { return m_handle->weak_count; }

      // non-std
      bool is_weak_ref() const { return m_is_weak_ref; }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
        reinterpret_cast<ElementType*>(m_handle->data), size())

      CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), size())

      void swap(shared_plain<ElementType>& other) {
        m_handle->swap(*(other.m_handle));
      }

      void reserve(const size_type& sz) {
        if (capacity() < sz) {
          shared_plain<ElementType> new_this(sz, reserve_flag());
          std::uninitialized_copy(begin(), end(), new_this.begin());
          new_this.m_set_size(size());
          new_this.swap(*this);
        }
      }

      // non-std
      shared_plain<ElementType>
      deep_copy() const {
        return shared_plain<ElementType>(begin(), end());
      }

      // non-std
      shared_plain<ElementType>
      weak_ref() const {
        return shared_plain<ElementType>(*this, weak_ref_flag());
      }

#     include <cctbx/array_family/push_back_etc.h>

    protected:

      size_type m_compute_new_capacity(const size_type& old_size,
                                       const size_type& n) {
        return old_size + std::max(old_size, n);
      }

      void m_insert_overflow(ElementType* pos,
                             const size_type& n, const ElementType& x,
                             bool at_end) {
        shared_plain<ElementType>
        new_this(m_compute_new_capacity(size(), n), reserve_flag());
        std::uninitialized_copy(begin(), pos, new_this.begin());
        new_this.m_set_size(pos - begin());
        if (n == 1) {
          new (new_this.end()) ElementType(x);
          new_this.m_incr_size(1);
        }
        else {
          std::uninitialized_fill_n(new_this.end(), n, x);
          new_this.m_incr_size(n);
        }
        if (!at_end) {
          std::uninitialized_copy(pos, end(), new_this.end());
          new_this.m_set_size(size() + n);
        }
        new_this.swap(*this);
      }

      void m_insert_overflow(ElementType* pos,
                             const ElementType* first,
                             const ElementType* last) {
        size_type n = last - first;
        shared_plain<ElementType>
        new_this(m_compute_new_capacity(size(), n), reserve_flag());
        std::uninitialized_copy(begin(), pos, new_this.begin());
        new_this.m_set_size(pos - begin());
        std::uninitialized_copy(first, last, new_this.end());
        new_this.m_incr_size(n);
        std::uninitialized_copy(pos, end(), new_this.end());
        new_this.m_set_size(size() + n);
        new_this.swap(*this);
      }

      void m_set_size(const size_type& sz) {
        m_handle->size = sz * element_size();
      }

      void m_incr_size(const size_type& n) {
        m_handle->size = (size() + n) * element_size();
      }

      void m_decr_size(const size_type& n) {
        m_handle->size = (size() - n) * element_size();
      }

      void m_dispose() {
        if (m_is_weak_ref) m_handle->weak_count--;
        else               m_handle->use_count--;
        if (m_handle->use_count == 0) {
          clear();
          if (m_handle->weak_count == 0) delete m_handle;
          else m_handle->deallocate();
        }
      }

      bool m_is_weak_ref;
      handle_type* m_handle;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_PLAIN_H
