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

#ifndef SCITBX_ARRAY_FAMILY_SHARED_PLAIN_H
#define SCITBX_ARRAY_FAMILY_SHARED_PLAIN_H

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/type_traits.h>

namespace scitbx { namespace af {

  struct weak_ref_flag {};

  namespace detail {
    const std::size_t global_max_size(static_cast<std::size_t>(-1));
  }


  class sharing_handle {
    public:
      sharing_handle()
        : use_count(1), weak_count(0), size(0), capacity(0),
          data(0)
      {}

      sharing_handle(weak_ref_flag)
        : use_count(0), weak_count(1), size(0), capacity(0),
          data(0)
      {}

      explicit
      sharing_handle(std::size_t const& sz)
        : use_count(1), weak_count(0), size(0), capacity(sz),
#if defined(SCITBX_ARRAY_FAMILY_SHARED_PLAIN_USE_STD_ALLOCATOR)
          data(alloc_.allocate(sz))
#else
          data(new char[sz])
#endif
      {}

      ~sharing_handle() {
#if defined(SCITBX_ARRAY_FAMILY_SHARED_PLAIN_USE_STD_ALLOCATOR)
        if (data) alloc_.deallocate(data, capacity);
#else
        delete[] data;
#endif
      }

      void deallocate() {
#if defined(SCITBX_ARRAY_FAMILY_SHARED_PLAIN_USE_STD_ALLOCATOR)
        alloc_.deallocate(data, capacity);
#else
        delete[] data;
#endif
        capacity = 0;
        data = 0;
      }

      void swap(sharing_handle& other) {
        std::swap(size, other.size);
        std::swap(capacity, other.capacity);
        std::swap(data, other.data);
      }

      std::size_t use_count;
      std::size_t weak_count;
      std::size_t size;
      std::size_t capacity;
      char* data;

    private:
#if defined(SCITBX_ARRAY_FAMILY_SHARED_PLAIN_USE_STD_ALLOCATOR)
      std::allocator<char> alloc_;
#endif
      sharing_handle(sharing_handle const&);
      sharing_handle& operator=(sharing_handle const&);
  };

  template <typename ElementType>
  class shared_plain
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      static size_type element_size() { return sizeof(ElementType); }

    public:
      shared_plain()
        : m_is_weak_ref(false),
          m_handle(new sharing_handle)
      {}

      explicit
      shared_plain(size_type const& sz)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle(sz * element_size()))
      {
        std::uninitialized_fill_n(begin(), sz, ElementType());
        m_handle->size = m_handle->capacity;
      }

      // non-std
      shared_plain(af::reserve const& sz)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle(sz() * element_size()))
      {}

      shared_plain(size_type const& sz, ElementType const& x)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle(sz * element_size()))
      {
        std::uninitialized_fill_n(begin(), sz, x);
        m_handle->size = m_handle->capacity;
      }

      // non-std
      template <typename FunctorType>
      shared_plain(size_type const& sz, init_functor<FunctorType> const& ftor)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle(sz * element_size()))
      {
        (*ftor.held)(begin(), sz);
        m_handle->size = m_handle->capacity;
      }

      shared_plain(const ElementType* first, const ElementType* last)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle((last - first) * element_size()))
      {
        std::uninitialized_copy(first, last, begin());
        m_handle->size = m_handle->capacity;
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      shared_plain(const OtherElementType* first, const OtherElementType* last)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle((last - first) * element_size()))
      {
        uninitialized_copy_typeconv(first, last, begin());
        m_handle->size = m_handle->capacity;
      }
#endif

      // non-std: shallow copy semantics
      shared_plain(shared_plain<ElementType> const& other)
        : m_is_weak_ref(other.m_is_weak_ref),
          m_handle(other.m_handle)
      {
        if (m_is_weak_ref) m_handle->weak_count++;
        else               m_handle->use_count++;
      }

      // non-std: shallow copy semantics, weak reference
      shared_plain(shared_plain<ElementType> const& other, weak_ref_flag)
        : m_is_weak_ref(true),
          m_handle(other.m_handle)
      {
        m_handle->weak_count++;
      }

      // non-std
      explicit
      shared_plain(sharing_handle* other_handle)
        : m_is_weak_ref(false),
          m_handle(other_handle)
      {
        SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->use_count++;
      }

      // non-std
      shared_plain(sharing_handle* other_handle, weak_ref_flag)
        : m_is_weak_ref(true),
          m_handle(other_handle)
      {
        SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->weak_count++;
      }

      // non-std
      template <typename OtherArrayType>
      shared_plain(array_adaptor<OtherArrayType> const& a_a)
        : m_is_weak_ref(false),
          m_handle(new sharing_handle)
      {
        OtherArrayType const& a = *(a_a.pointee);
        reserve(a.size());
        for(std::size_t i=0;i<a.size();i++) push_back(a[i]);
      }

      ~shared_plain() {
        m_dispose();
      }

      // non-std: shallow copy semantics
      shared_plain<ElementType>&
      operator=(shared_plain<ElementType> const& other)
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

      // non-std
      std::size_t id() const
      {
        return reinterpret_cast<std::size_t>(m_handle);
      }

      // non-std
      sharing_handle* handle() {
        SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
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

      SCITBX_ARRAY_FAMILY_BEGIN_END_ETC(shared_plain,
        reinterpret_cast<ElementType*>(m_handle->data), size())

      SCITBX_ARRAY_FAMILY_TAKE_REF(begin(), size())

      void swap(shared_plain<ElementType>& other) {
        m_handle->swap(*(other.m_handle));
      }

      void reserve(size_type const& sz) {
        if (capacity() < sz) {
          shared_plain<ElementType> new_this(((af::reserve(sz))));
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

#     include <scitbx/array_family/detail/push_back_etc.h>

    protected:

      size_type m_compute_new_capacity(size_type const& old_size,
                                       size_type const& n) {
        return old_size + std::max(old_size, n);
      }

      void m_insert_overflow(ElementType* pos,
                             size_type const& n, ElementType const& x,
                             bool at_end) {
        shared_plain<ElementType>
          new_this((af::reserve(m_compute_new_capacity(size(), n))));
        std::uninitialized_copy(begin(), pos, new_this.begin());
        new_this.m_set_size(pos - begin());
        if (n == 1) {
          new (new_this.end()) ElementType(x);
          new_this.m_incr_size(1);
        }
        else {
          // The next call is known to lead to a segmentation fault under
          // RedHat 8.0 if there is not enough memory (observed with
          // ElementType = std::vector<vector_element> and
          // struct vector_element { unsigned i; unsigned j; };).
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
          new_this((af::reserve(m_compute_new_capacity(size(), n))));
        std::uninitialized_copy(begin(), pos, new_this.begin());
        new_this.m_set_size(pos - begin());
        std::uninitialized_copy(first, last, new_this.end());
        new_this.m_incr_size(n);
        std::uninitialized_copy(pos, end(), new_this.end());
        new_this.m_set_size(size() + n);
        new_this.swap(*this);
      }

      void m_set_size(size_type const& sz) {
        m_handle->size = sz * element_size();
      }

      void m_incr_size(size_type const& n) {
        m_handle->size = (size() + n) * element_size();
      }

      void m_decr_size(size_type const& n) {
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
      sharing_handle* m_handle;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SHARED_PLAIN_H
