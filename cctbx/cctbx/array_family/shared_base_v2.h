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

#ifndef CCTBX_ARRAY_FAMILY_SHARED_BASE_V2_H
#define CCTBX_ARRAY_FAMILY_SHARED_BASE_V2_H

#include <algorithm>
#include <cctbx/array_family/ref.h>

namespace stlp {

  template <class T1, class T2>
  inline void construct(T1* p, const T2& val) {
    new (p) T1(val);
  }

  template <class _Value, class _Tp>
  class alloc_proxy
  public:
    _Value m_data;

    inline _Tp* allocate(size_t __n) { 
      return 0;
    }
    inline void deallocate(_Tp* __p, size_t __n) { 
    }
  };

}

namespace cctbx { namespace af {
  namespace detail {

    template <typename ElementType,
              std::size_t SmallestAutoCapacity = 8>
    class basic_storage {
      public:
        typedef std::size_t size_type;

        basic_storage()
          : m_capacity(0), m_size(0), m_use_count(1),
            m_data(0)
        {}

        explicit
        basic_storage(const size_type& sz)
          : m_capacity(sz), m_size(sz), m_use_count(1),
            m_data(new ElementType[sz])
        {}

        ~basic_storage() {
          delete[] m_data;
        }

        size_type size() const { return m_size; }
        size_type capacity() const { return m_capacity; }

              long& use_count()       { return m_use_count; }
        const long& use_count() const { return m_use_count; }

        // no const: data are not considered part of the type
        ElementType* begin() const { return m_data; }

        void swap(basic_storage<ElementType>& other){
          std::swap(*this, other);
        }

        void reserve(const size_type& new_capacity) {
          if (new_capacity > m_capacity) {
            m_capacity = new_capacity;
            ElementType* new_data = new ElementType[m_capacity];
            std::copy(m_data, m_data + m_size, new_data);
            delete [] m_data;
            m_data = new_data;
          }
        }

        void resize(const size_type& new_size) {
          reserve(new_size);
          m_size = new_size;
        }

        void auto_resize(const size_type& new_size) {
          if (new_size > 0) {
            reserve(std::max(std::max(
              new_size, SmallestAutoCapacity), m_capacity * 2));
          }
          m_size = new_size;
        }

      private:
        size_type m_capacity;
        size_type m_size;
        long m_use_count;
        ElementType* m_data;
    };

  } // namespace detail

  template <class ElementType>
  class shared_base_v2
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

    protected:
      struct m_false_type {};
      struct m_true_type {};

      typedef typename m_false_type m_has_trivial_assignment_operator;
      typedef typename m_false_type m_is_pod_type;

    public:
            ElementType* begin()       { return this->m_start; }
      const ElementType* begin() const { return this->m_start; }
            ElementType* end()        { return this->m_finish; }
      const ElementType* end() const  { return this->m_finish; }

      size_type size() const {
        return size_type(this->m_finish - this->m_start);
      }

      size_type max_size() const {
        return size_type(-1) / sizeof(ElementType);
      }

      size_type capacity() const {
        return size_type(this->m_end_of_storage.m_data - this->m_start);
      }

      bool empty() const {
        return this->m_start == this->m_finish;
      }

            ElementType& operator[](size_type i) {
        return *(begin() + i);
      }
      const ElementType& operator[](size_type i) const {
        return *(begin() + i);
      }

            ElementType& front()       { return *begin(); }
      const ElementType& front() const { return *begin(); }
            ElementType& back()       { return *(end() - 1); }
      const ElementType& back() const { return *(end() - 1); }

            ElementType& at(const size_type& i) {
        if (i >= sz) throw_range_error();
        return (*this)[i];
      }
      const ElementType& at(const size_type& i) const {
        if (i >= sz) throw_range_error();
        return (*this)[i];
      }

      shared_base_v2()
        : m_ptr_basic_storage(new basic_storage)
      {}

      shared_base_v2(size_type sz, const ElementType& x)
        : m_ptr_basic_storage(new basic_storage(sz * sizeof(ElementType)))
      {
        this->m_finish = std::uninitialized_fill_n(this->m_start, sz, x);
      }

      explicit
      shared_base_v2(size_type sz)
        : m_ptr_basic_storage(new basic_storage(sz * sizeof(ElementType)))
      {
        this->m_finish = std::uninitialized_fill_n(
          this->m_start, sz, ElementType());
      }

      shared_base_v2(const shared_base_v2<ElementType>& other)
        : m_ptr_basic_storage(
            new basic_storage(other.size() * sizeof(ElementType)))
      {
        this->m_finish = std::uninitialized_copy(
          other.begin(), other.end(), this->m_start);
      }

      shared_base_v2(const ElementType* first, const ElementType* last)
        : m_ptr_basic_storage(
            new basic_storage((last - first) * sizeof(ElementType)))
      {
        size_type sz = last - first;
        this->m_start = this->m_end_of_storage.allocate(sz);
        this->m_end_of_storage.m_data = this->m_start + sz;
        this->m_finish = std::uninitialized_copy(
          first, last, this->m_start);
      }

      ~shared_base_v2() {
        detail::destroy_array_elements(begin(), end());
      }

      shared_base_v2<ElementType>&
      operator=(const shared_base_v2<ElementType>& other)
      {
        if (&other != this) {
          const size_type other_size = other.size();
          if (other_size > capacity()) {
            ElementType* tmp = m_allocate_and_copy(other_size,
              (const ElementType*)other.m_start+0,
              (const ElementType*)other.m_finish+0);
            m_clear();
            this->m_start = tmp;
            this->m_end_of_storage.m_data = this->m_start + other_size;
          }
          else if (size() >= other_size) {
            ElementType* i = std::copy(
              other.m_start, other.m_finish, this->m_start);
            detail::destroy_array_elements(i, this->m_finish);
          }
          else {
            std::copy(
              other.m_start, other.m_start + size(), this->m_start);
            std::uninitialized_copy(
              other.m_start + size(), other.m_finish, this->m_finish);
          }
          this->m_finish = this->m_start + other_size;
        }
        return *this;
      }

      void swap(shared_base_v2<ElementType>& other) {
        std::swap(this->m_start, other.m_start);
        std::swap(this->m_finish, other.m_finish);
        std::swap(this->m_end_of_storage, other.m_end_of_storage);
      }

      reserve(size_type sz) {
        if (capacity() < sz) {
          const size_type old_size = size();
          ElementType* tmp;
          if (this->m_start) {
            tmp = m_allocate_and_copy(sz, this->m_start, this->m_finish);
            m_clear();
          }
          else {
            tmp = this->m_end_of_storage.allocate(sz);
          }
          m_set(tmp, tmp + old_size, tmp + sz);
        }
      }

      void assign(size_type sz, const ElementType& x) {
        m_fill_assign(sz, x);
      }

      void assign(const ElementType* first, const ElementType* last)
      {
        size_type n = last - first;
        if (n > capacity()) {
          ElementType* tmp = m_allocate_and_copy(n, first, last);
          m_clear();
          m_set(tmp, tmp + n, tmp + n);
        }
        else if (size() >= n) {
          ElementType* new_finish = std::copy(
            first, last, this->m_start);
          detail::destroy_array_elements(new_finish, this->m_finish);
          this->m_finish = new_finish;
        }
        else {
          const ElementType* mid = first + size() ;
          std::copy(first, mid, this->m_start);
          this->m_finish = std::uninitialized_copy(
            mid, last, this->m_finish);
        }
      }

      void push_back(const ElementType& x) {
        if (this->m_finish != this->m_end_of_storage.m_data) {
          stlp::construct(this->m_finish, x);
          ++this->m_finish;
        }
        else {
          m_insert_overflow(this->m_finish, x, m_is_pod_type(), 1UL, true);
        }
      }

      ElementType* insert(ElementType* position, const ElementType& x) {
        size_type n = position - begin();
        if (this->m_finish != this->m_end_of_storage.m_data) {
          if (position == end()) {
            stlp::construct(this->m_finish, x);
            ++this->m_finish;
          }
          else {
            stlp::construct(this->m_finish, *(this->m_finish - 1));
            ++this->m_finish;
            ElementType x_copy = x;
            std::copy_backward(
              position, this->m_finish - 2, this->m_finish - 1);
            *position = x_copy;
          }
        }
        else {
          m_insert_overflow(position, x, m_is_pod_type(), 1UL);
        }
        return begin() + n;
      }

      void insert(ElementType* position,
                  const ElementType* first, const ElementType* last)
      {
        if (first != last) {
          size_type n = last - first;
          if (this->m_end_of_storage.m_data - this->m_finish >= n) {
            const size_type elems_after = this->m_finish - position;
            ElementType* old_finish = this->m_finish;
            if (elems_after > n) {
              std::uninitialized_copy(
                this->m_finish - n, this->m_finish, this->m_finish);
              this->m_finish += n;
              std::copy_backward(position, old_finish - n, old_finish);
              std::copy(first, last, position);
            }
            else {
              const ElementType* mid = first + elems_after;
              std::uninitialized_copy(
                mid, last, this->m_finish);
              this->m_finish += n - elems_after;
              std::uninitialized_copy(
                position, old_finish, this->m_finish);
              this->m_finish += elems_after;
              std::copy(first, mid, position);
            }
          }
          else {
            const size_type old_size = size();
            const size_type len = old_size + std::max(old_size, n);
            ElementType* new_start = this->m_end_of_storage.allocate(len);
            ElementType* new_finish = new_start;
            try {
              new_finish = std::uninitialized_copy(
                this->m_start, position, new_start);
              new_finish = std::uninitialized_copy(
                first, last, new_finish);
              new_finish = std::uninitialized_copy(
                position, this->m_finish, new_finish);
            }
            catch (...) {
              detail::destroy_array_elements(new_start, new_finish);
              this->m_end_of_storage.deallocate(new_start, len);
              throw;
            }
            m_clear();
            m_set(new_start, new_finish, new_start + len);
          }
        }
      }

      void insert(ElementType* pos, size_type n, const ElementType& x) {
        m_fill_insert(pos, n, x);
      }

      void pop_back() {
        --this->m_finish;
        detail::destroy_array_element(this->m_finish);
      }

      ElementType* erase(ElementType* position) {
        if (position + 1 != end()) {
          std::copy(position + 1, this->m_finish, position);
        }
        --this->m_finish;
        detail::destroy_array_element(this->m_finish);
        return position;
      }

      ElementType* erase(ElementType* first, ElementType* last) {
        ElementType* i = std::copy(last, this->m_finish, first);
        detail::destroy_array_elements(i, this->m_finish);
        this->m_finish = i;
        return first;
      }

      void resize(size_type new_size, const ElementType& x) {
        if (new_size < size())  {
          erase(begin() + new_size, end());
        }
        else {
          insert(end(), new_size - size(), x);
        }
      }

      void resize(size_type new_size) {
        resize(new_size, ElementType());
      }

      void clear() {
        erase(begin(), end());
      }

    protected:

      // for non-POD types
      void m_insert_overflow(ElementType* position, const ElementType& x,
                             const m_false_type&,
                             size_type fill_len, bool at_end = false) {
        const size_type old_size = size();
        const size_type len = old_size + std::max(old_size, fill_len);

        ElementType* new_start = this->m_end_of_storage.allocate(len);
        ElementType* new_finish = new_start;
        try {
          new_finish = std::uninitialized_copy(
            this->m_start, position, new_start);
          // handle insertion
          if (fill_len == 1) {
            stlp::construct(new_finish, x);
            ++new_finish;
          }
          else {
            new_finish = std::uninitialized_fill_n(
              new_finish, fill_len, x);
          }
          if (!at_end)
            // copy remainder
            new_finish = std::uninitialized_copy(
              position, this->m_finish, new_finish);
        }
        catch (...) {
          detail::destroy_array_elements(new_start, new_finish);
          this->m_end_of_storage.deallocate(new_start, len);
          throw;
        }
        m_clear();
        m_set(new_start, new_finish, new_start + len);
      }

      // for POD types
      void m_insert_overflow(ElementType* position, const ElementType& x,
                             const m_true_type&,
                             size_type fill_len, bool at_end = false) {
        const size_type old_size = size();
        const size_type len = old_size + std::max(old_size, fill_len);

        ElementType* new_start = this->m_end_of_storage.allocate(len);
        ElementType* new_finish = (ElementType*) std::copy(
          this->m_start, position, new_start);
          // handle insertion
        new_finish = fill_n(new_finish, fill_len, x);
        if (!at_end)
          // copy remainder
          new_finish = (ElementType*) std::copy(
            position, this->m_finish, new_finish);
        m_clear();
        m_set(new_start, new_finish, new_start + len);
      }

      void m_clear() {
        detail::destroy_array_elements(this->m_start, this->m_finish);
        this->m_end_of_storage.deallocate(
          this->m_start, this->m_end_of_storage.m_data - this->m_start);
      }

      void m_set(ElementType* s, ElementType* f, ElementType* e) {
        this->m_start = s;
        this->m_finish = f;
        this->m_end_of_storage.m_data = e;
      }

      ElementType* m_allocate_and_copy(
        size_type n, const ElementType* first, const ElementType* last)
      {
        ElementType* result = this->m_end_of_storage.allocate(n);
        try {
          std::uninitialized_copy(first, last, result);
          return result;
        }
        catch (...) {
          this->m_end_of_storage.deallocate(result, n);
          throw;
        }
      }

      m_fill_insert(ElementType* position, size_type n, const ElementType& x) {
        if (n != 0) {
          if (size_type(this->m_end_of_storage.m_data - this->m_finish) >= n) {
            ElementType x_copy = x;
            const size_type elems_after = this->m_finish - position;
            ElementType* old_finish = this->m_finish;
            if (elems_after > n) {
              std::uninitialized_copy(
                this->m_finish - n, this->m_finish, this->m_finish);
              this->m_finish += n;
              std::copy_backward(position, old_finish - n, old_finish);
              std::fill(position, position + n, x_copy);
            }
            else {
              std::uninitialized_fill_n(
                this->m_finish, n - elems_after, x_copy);
              this->m_finish += n - elems_after;
              std::uninitialized_copy(position, old_finish, this->m_finish);
              this->m_finish += elems_after;
              std::fill(position, old_finish, x_copy);
            }
          }
          else {
            m_insert_overflow(position, x, m_is_pod_type(), n);
          }
        }
      }

      m_fill_assign(size_t n, const ElementType& val) {
        if (n > capacity()) {
          shared_base_v2<ElementType> tmp(n, val, get_allocator());
          tmp.swap(*this);
        }
        else if (n > size()) {
          fill(begin(), end(), val);
          this->m_finish = std::uninitialized_fill_n(
            this->m_finish, n - size(), val);
        }
        else {
          erase(std::fill_n(begin(), n, val), end());
        }
      }

      detail::basic_storage* m_ptr_basic_storage;
      ElementType* m_start;
      ElementType* m_finish;
      stlp::alloc_proxy<_Tp*, _Tp> m_end_of_storage;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_SHARED_BASE_V2_H
