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

#include <scitbx/error.h>

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/type_traits.h>

#include "polymorphic_allocator.h"

namespace scitbx { namespace af {

  struct weak_ref_flag {};

  namespace detail {
    const std::size_t global_max_size(static_cast<std::size_t>(-1));
  }

  /** Data store that is shared between shared_plain objects.
   *
   * Size is allocated in bytes, and it contains it's own reference counts.
   * Although size and capacity is determined on creation, the inner
   * configuration of this class is managed by the shared_plain owners.
   */
  class sharing_handle {
    public:
      // General allocator - should be std::byte once C++11 is supported
      typedef pmr::polymorphic_allocator<unsigned char> allocator_type;

      sharing_handle(allocator_type a = allocator_type())
        : use_count(1), weak_count(0), size(0), capacity(0),
          data(NULL), m_allocator(a)
      {}

      sharing_handle(weak_ref_flag, allocator_type a = allocator_type())
        : use_count(0), weak_count(1), size(0), capacity(0),
          data(NULL), m_allocator(a)
      {}

      explicit
      sharing_handle(std::size_t const& sz, allocator_type a = allocator_type())
        : use_count(1), weak_count(0), size(0), capacity(sz),
          data(NULL), m_allocator(a)
      {
        // Allow calling with size zero to match other constructors
        if (sz != 0) data = m_allocator.allocate(sz);
      }

      ~sharing_handle() {
        if (data) m_allocator.deallocate(data, capacity);
        data = NULL;
      }

      void deallocate() {
        if (data) m_allocator.deallocate(data, capacity);
        capacity = 0;
        data = NULL;
      }

      void swap(sharing_handle& other) {
        // Swapping doesn't make sense if the allocators don't match -
        // it would be a copy instead, in which case we almost definitely
        // won't want a symmetric copy.
        SCITBX_ASSERT(m_allocator == other.m_allocator);

        std::swap(size, other.size);
        std::swap(capacity, other.capacity);
        std::swap(data, other.data);
      }

      std::size_t use_count;
      std::size_t weak_count;
      std::size_t size;
      std::size_t capacity;
      allocator_type::value_type* data;

      allocator_type get_allocator() const { return m_allocator; }

    private:
      // No copy or assignment operators
      sharing_handle(sharing_handle const&);
      sharing_handle& operator=(sharing_handle const&);

      allocator_type m_allocator;
  };

  /**
   * Container that shares data between it's instances.
   *
   * Conceptually, acts like a shared_ptr<vector<ElementType>>, but
   * without the need to use the dereferencing operator to operate
   * on the vector.
   *
   * Has additional functionality to create weak references to the
   * main data set.
   */
  template <typename ElementType>
  class shared_plain
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef pmr::polymorphic_allocator<unsigned char> allocator_type;

      static size_type element_size() { return sizeof(ElementType); }

    public:
      explicit
      shared_plain(allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(0, m_allocator))
      {}

      explicit
      shared_plain(size_type const& sz, allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(sz, m_allocator))
      {
        std::uninitialized_fill_n(begin(), sz, ElementType());
        m_handle->size = m_handle->capacity;
      }

      // non-std
      shared_plain(af::reserve const& sz, allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(sz(), m_allocator))
      { }

      shared_plain(size_type const& sz, ElementType const& x, allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(sz, m_allocator))
      {
        std::uninitialized_fill_n(begin(), sz, x);
        m_handle->size = m_handle->capacity;
      }

      // non-std
      template <typename FunctorType>
      shared_plain(size_type const& sz, init_functor<FunctorType> const& ftor, allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(sz, m_allocator))
      {
        (*ftor.held)(begin(), sz);
        m_handle->size = m_handle->capacity;
      }

      /// Create a shared_plain by copying out of an existing memory buffer
      shared_plain(const ElementType* first, const ElementType* last, allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(last-first, m_allocator))
      {
        std::uninitialized_copy(first, last, begin());
        m_handle->size = m_handle->capacity;
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      template <typename OtherElementType>
      shared_plain(const OtherElementType* first, const OtherElementType* last, allocator_type a = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(a),
          m_handle(m_allocate_handle(last-first, m_allocator))
      {
        uninitialized_copy_typeconv(first, last, begin());
        m_handle->size = m_handle->capacity;
      }
#endif

      // non-std: shallow copy semantics
      // Allocator: Normally copy would give a chance to change allocator -
      //            but we're way out of standard behaviour territory here,
      //            so propogate the other allocator to this instance so that
      //            if we ever need to do a reallocate (which at time of
      //            writing isn't done), then it's on the correct resource
      shared_plain(shared_plain<ElementType> const& other)
        : m_is_weak_ref(other.m_is_weak_ref),
          m_allocator(other.m_allocator),
          m_handle(other.m_handle)
      {
        if (m_is_weak_ref) m_handle->weak_count++;
        else               m_handle->use_count++;
      }

      // non-std: shallow copy semantics, weak reference
      shared_plain(shared_plain<ElementType> const& other, weak_ref_flag)
        : m_is_weak_ref(true),
          m_allocator(other.m_allocator),
          m_handle(other.m_handle)
      {
        m_handle->weak_count++;
      }

      // non-std
      // Allocator: Tricky here, as technically we don't know what allocator
      //            the other handle was allocated on, we can't deallocate
      //            it reliably. However, in all current functionality
      //            it's a safe bet that it was allocated on the same
      //            allocator that the handle is using, so use that and
      //            hope we don't cause any segmentation faults/leaks.
      explicit
      shared_plain(sharing_handle* other_handle)
          : m_is_weak_ref(false),
            m_allocator(other_handle->get_allocator()),
            m_handle(other_handle) {
        SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->use_count++;
      }

      // non-std
      // Allocator issue same as non-weak version. Not advised.
      shared_plain(sharing_handle* other_handle, weak_ref_flag)
          : m_is_weak_ref(true),
            m_allocator(other_handle->get_allocator()),
            m_handle(other_handle) {
        SCITBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        m_handle->weak_count++;
      }

      // non-std
      template <typename OtherArrayType>
      shared_plain(array_adaptor<OtherArrayType> const& a_a, allocator_type allocator = allocator_type())
        : m_is_weak_ref(false),
          m_allocator(allocator),
          m_handle(m_allocate_handle(0, m_allocator))
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
          // Non-copy assign only works without copy for identical allocators
          SCITBX_ASSERT(m_allocator == other.m_allocator);

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
          shared_plain<ElementType> new_this(((af::reserve(sz))), m_allocator);
          std::uninitialized_copy(begin(), end(), new_this.begin());
          new_this.m_set_size(size());
          new_this.swap(*this);
        }
      }

      // non-std
      shared_plain<ElementType>
      deep_copy(allocator_type a) const {
        return shared_plain<ElementType>(begin(), end(), a);
      }

      shared_plain<ElementType>
      deep_copy() const {
        return deep_copy(allocator_type());
      }

      // non-std
      shared_plain<ElementType>
      weak_ref() const {
        return shared_plain<ElementType>(*this, weak_ref_flag());
      }

#     include <scitbx/array_family/detail/push_back_etc.h>

      allocator_type get_allocator() const { return m_allocator; }

    protected:

      size_type m_compute_new_capacity(size_type const& old_size,
                                       size_type const& n) {
        return old_size + std::max(old_size, n);
      }

      void m_insert_overflow(ElementType* pos,
                             size_type const& n, ElementType const& x,
                             bool at_end) {
        shared_plain<ElementType>
          new_this((af::reserve(m_compute_new_capacity(size(), n))), m_allocator);
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
          new_this((af::reserve(m_compute_new_capacity(size(), n))), m_allocator);
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
          if (m_handle->weak_count == 0) {
            // Completely delete the shared handle now
            m_allocator.destroy(m_handle);
            m_allocator.deallocate((allocator_type::pointer)m_handle, sizeof(sharing_handle));
            m_handle = 0;
          }
          else {
            // Keep the handle around, but delete it's internal memory
            m_handle->deallocate();
          }
        }
      }

      bool m_is_weak_ref;
      allocator_type m_allocator;
      sharing_handle* m_handle;

    private:
      /** Convenience function to allocate and create a sharing handle
       *
       *  Because we have many, many constructors this would become overly
       *  repetitious. With C++11 this would probably be better fulfilled
       *  by delegating constructors.
       */
      static sharing_handle* m_allocate_handle(size_t size,
                                               allocator_type& allocator) {
        sharing_handle* handle =
            (sharing_handle*)allocator.allocate(sizeof(sharing_handle));
        allocator.construct(handle, size * element_size(), allocator);
        return handle;
     }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_SHARED_PLAIN_H
