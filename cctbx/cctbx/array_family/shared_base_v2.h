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

// XXX weak reference count

#ifndef CCTBX_ARRAY_FAMILY_SHARED_BASE_V2_H
#define CCTBX_ARRAY_FAMILY_SHARED_BASE_V2_H

#include <algorithm>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/type_traits.h>

namespace cctbx { namespace af {

  namespace detail {
    const std::size_t global_max_size(std::size_t(-1));
  }

  template <typename ElementType>
  class shared_base_v2
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      static size_type element_size() { return sizeof(ElementType); }

    protected:
      class handle_type {
        public:
          handle_type()
            : use_count(1), size(0), capacity(0),
              data(0)
          {}

          explicit
          handle_type(const size_type& sz)
            : use_count(1), size(0), capacity(sz),
              data(new char[sz])
          {}

          ~handle_type() {
            delete[] data;
          }

          void swap(handle_type& other) {
            std::swap(size, other.size);
            std::swap(capacity, other.capacity);
            std::swap(data, other.data);
          }

          size_type use_count;
          size_type size;
          size_type capacity;
          char* data;

        private:
          handle_type(const handle_type&);
          handle_type& operator=(const handle_type&);
      };

    public:
      shared_base_v2()
        : m_handle(new handle_type)
      {}

      explicit
      shared_base_v2(const size_type& sz)
        : m_handle(new handle_type(sz * element_size()))
      {
        std::uninitialized_fill_n(begin(), sz, ElementType());
        m_handle->size = m_handle->capacity;
      }

      // non-std
      shared_base_v2(const size_type& sz, reserve_flag)
        : m_handle(new handle_type(sz * element_size()))
      {}

      shared_base_v2(const size_type& sz, const ElementType& x)
        : m_handle(new handle_type(sz * element_size()))
      {
        std::uninitialized_fill_n(begin(), sz, x);
        m_handle->size = m_handle->capacity;
      }

      shared_base_v2(const ElementType* first, const ElementType* last)
        : m_handle(new handle_type((last - first) * element_size()))
      {
        std::uninitialized_copy(first, last, begin());
        m_handle->size = m_handle->capacity;
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
                                        // 1200 == VC++ 6.0
      template <typename OtherElementType>
      shared_base_v2(
        const OtherElementType* first, const OtherElementType* last)
        : m_handle(new handle_type((last - first) * element_size()))
      {
        uninitialized_copy_typeconv(first, last, begin());
        m_handle->size = m_handle->capacity;
      }
#endif

      // non-std: shallow copy semantics
      shared_base_v2(const shared_base_v2<ElementType>& other)
        : m_handle(other.m_handle)
      {
        m_handle->use_count++;
      }

      explicit
      shared_base_v2(const handle_type& handle)
        : m_handle(handle)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
      }

      ~shared_base_v2() {
        m_dispose();
      }

      // non-std: shallow copy semantics
      shared_base_v2<ElementType>&
      operator=(const shared_base_v2<ElementType>& other)
      {
        if (m_handle != other.m_handle) {
          m_dispose();
          m_handle = other.m_handle;
          m_handle->use_count++;
        }
        return *this;
      }

      // non-std
            handle_type& handle()       {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        return m_handle;
      }
      // non-std
      const handle_type& handle() const {
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

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
        reinterpret_cast<ElementType*>(m_handle->data), size())

      CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), size())

      void swap(shared_base_v2<ElementType>& other) {
        m_handle->swap(*(other.m_handle));
      }

      void reserve(const size_type& sz) {
        if (capacity() < sz) {
          shared_base_v2<ElementType> new_this(sz, reserve_flag());
          new_this.assign(begin(), end());
          new_this.swap(*this);
        }
      }

      void assign(const size_type& sz, const ElementType& x) {
        if (sz > capacity()) {
          shared_base_v2<ElementType> new_this(sz, x);
          new_this.swap(*this);
        }
        else if (sz > size()) {
          std::fill(begin(), end(), x);
          std::uninitialized_fill(end(), begin() + sz, x);
          m_set_size(sz);
        }
        else {
          std::fill_n(begin(), sz, x);
          erase(begin() + sz, end());
        }
      }

      void assign(const ElementType* first, const ElementType* last)
      {
        size_type sz = last - first;
        if (sz > capacity()) {
          shared_base_v2<ElementType> new_this(first, last);
          new_this.swap(*this);
        }
        else if (sz > size()) {
          const ElementType* mid = first + size() ;
          std::copy(first, mid, begin());
          std::uninitialized_copy(mid, last, end());
          m_set_size(sz);
        }
        else {
          std::copy(first, last, begin());
          erase(begin() + sz, end());
        }
      }

      void push_back(const ElementType& x) {
        if (size() < capacity()) {
          new (end()) ElementType(x);
          m_incr_size(1);
        }
        else {
          m_insert_overflow(end(), x, size_type(1), true);
        }
      }

      void pop_back() {
        m_decr_size(1);
        detail::destroy_array_element(end());
      }

      ElementType* insert(ElementType* pos, const ElementType& x) {
        size_type n = pos - begin();
        if (size() == capacity()) {
          m_insert_overflow(pos, x, size_type(1), false);
        }
        else {
          if (pos == end()) {
            new (end()) ElementType(x);
            m_handle->size = (size() + 1) * element_size();
          }
          else {
            new (end()) ElementType(*(end() - 1));
            m_incr_size(1);
            ElementType x_copy = x;
            std::copy_backward(pos, end() - 2, end() - 1);
            *pos = x_copy;
          }
        }
        return begin() + n;
      }

      // restricted to ElementType
      void insert(ElementType* pos,
                  const ElementType* first, const ElementType* last)
      {
        size_type n = last - first;
        if (n == 0) return;
        size_type new_size = size() + n;
        if (new_size <= capacity()) {
          size_type n_move_up = end() - pos;
          ElementType* old_end = end();
          if (n_move_up > n) {
            std::uninitialized_copy(end() - n, end(), end());
            m_incr_size(n);
            std::copy_backward(pos, old_end - n, old_end);
            std::copy(first, last, pos);
          }
          else {
            const ElementType* mid = first + n_move_up;
            std::uninitialized_copy(mid, last, end());
            m_incr_size(n - n_move_up);
            std::uninitialized_copy(pos, old_end, end());
            m_incr_size(n_move_up);
            std::copy(first, mid, pos);
          }
        }
        else {
          size_type old_size = size();
          size_type new_capacity = old_size + std::max(old_size, n);
          shared_base_v2<ElementType>
          new_this(new_capacity, reserve_flag());
          std::uninitialized_copy(begin(), pos, new_this.begin());
          new_this.m_set_size(pos - begin());
          std::uninitialized_copy(first, last, new_this.end());
          new_this.m_incr_size(n);
          std::uninitialized_copy(pos, end(), new_this.end());
          new_this.m_set_size(new_size);
          new_this.swap(*this);
        }
      }

      void insert(ElementType* pos, const size_type& n, const ElementType& x) {
        if (n == 0) return;
        if (size() + n > capacity()) {
          m_insert_overflow(pos, x, n, false);
        }
        else {
          ElementType x_copy = x;
          size_type n_move_up = end() - pos;
          ElementType* old_end = end();
          if (n_move_up > n) {
            std::uninitialized_copy(end() - n, end(), end());
            m_incr_size(n);
            std::copy_backward(pos, old_end - n, old_end);
            std::fill_n(pos, n, x_copy);
          }
          else {
            std::uninitialized_fill_n(end(), n - n_move_up, x_copy);
            m_incr_size(n - n_move_up);
            std::uninitialized_copy(pos, old_end, end());
            m_incr_size(n_move_up);
            std::fill(pos, old_end, x_copy);
          }
        }
      }

      ElementType* erase(ElementType* pos) {
        if (pos + 1 != end()) {
          std::copy(pos + 1, end(), pos);
        }
        m_decr_size(1);
        detail::destroy_array_element(end());
        return pos;
      }

      ElementType* erase(ElementType* first, ElementType* last) {
        ElementType* i = std::copy(last, end(), first);
        detail::destroy_array_elements(i, end());
        m_decr_size(last - first);
        return first;
      }

      void resize(const size_type& new_size, const ElementType& x) {
        if (new_size < size())  {
          erase(begin() + new_size, end());
        }
        else {
          insert(end(), new_size - size(), x);
        }
      }

      void resize(const size_type& new_size) {
        resize(new_size, ElementType());
      }

      void clear() {
        erase(begin(), end());
      }

      // non-std
      shared_base_v2<ElementType>
      deep_copy() const {
        shared_base_v2<ElementType> result(size(), reserve_flag());
        result.assign(begin(), end());
        return result;
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
                                        // 1200 == VC++ 6.0
      // non-std
      template <typename OtherElementType>
      void
      assign(const OtherElementType* first, const OtherElementType* last) {
        clear();
        size_type sz = last - first;
        reserve(sz);
        uninitialized_copy_typeconv(first, last, begin());
        m_set_size(sz);
      }
#endif

      // non-std
      template <typename OtherArrayType>
      void
      assign(const OtherArrayType& other) {
        assign(other.begin(), other.end());
      }

    protected:

      void m_insert_overflow(ElementType* pos, const ElementType& x,
                             const size_type& n, bool at_end) {
        size_type old_size = size();
        size_type new_capacity = old_size + std::max(old_size, n);
        shared_base_v2<ElementType>
        new_this(new_capacity, reserve_flag());
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
          new_this.m_set_size(old_size + n);
        }
        new_this.swap(*this);
      }

      void m_dispose() {
        if (m_handle->use_count > 1) {
          m_handle->use_count--;
        }
        else {
          erase(begin(), end());
          delete m_handle;
        }
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

      handle_type* m_handle;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_BASE_V2_H
