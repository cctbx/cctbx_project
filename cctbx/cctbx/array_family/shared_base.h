// $Id$
/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Reorganized (rwgk)
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_BASE_H
#define CCTBX_ARRAY_FAMILY_SHARED_BASE_H

#include <algorithm>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/type_traits.h>

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

    template <typename ElementType>
    class managed_storage {
      public:
        typedef ElementType element_type;
        typedef std::size_t size_type;
        typedef basic_storage<element_type> basic_storage_type;

        explicit
        managed_storage(const size_type& sz)
          : m_ptr_basic_storage(new basic_storage_type(sz))
        {}

        managed_storage(const managed_storage& other)
          : m_ptr_basic_storage(other.m_ptr_basic_storage) {
          m_ptr_basic_storage->use_count()++;
        }

        managed_storage& operator=(const managed_storage& other) {
          if (m_ptr_basic_storage != other.m_ptr_basic_storage) {
            dispose();
            m_ptr_basic_storage = other.m_ptr_basic_storage;
            m_ptr_basic_storage->use_count()++;
          }
          return *this;
        }

        ~managed_storage() { dispose(); }

        void swap(managed_storage<element_type>& other) {
          std::swap(m_ptr_basic_storage, other.m_ptr_basic_storage);
          std::swap(      m_ptr_basic_storage->use_count(),
                    other.m_ptr_basic_storage->use_count());
        }

        size_type size() const { return m_ptr_basic_storage->size(); }
        size_type capacity() const { return m_ptr_basic_storage->capacity(); }
        long use_count() const { return m_ptr_basic_storage->use_count(); }

        // no const: data are not considered part of the type
        element_type* begin() const { return m_ptr_basic_storage->begin(); }

        void reserve(const size_type& new_size) {
          m_ptr_basic_storage->reserve(new_size);
        }

        void resize(const size_type& new_size) {
          m_ptr_basic_storage->resize(new_size);
        }

        void auto_resize(const size_type& new_size) {
          m_ptr_basic_storage->auto_resize(new_size);
        }

      private:
        basic_storage_type* m_ptr_basic_storage;

        void dispose() {
          m_ptr_basic_storage->use_count()--;
          if (m_ptr_basic_storage->use_count() == 0) {
            delete m_ptr_basic_storage;
          }
        }
    };

  } //namespace detail

  template <typename ElementType>
  class shared_base
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef detail::managed_storage<char> handle_type;

      static size_type element_size() { return sizeof(ElementType); }

      explicit
      shared_base(const size_type& sz = 0,
                  const ElementType& x = ElementType())
        : m_handle(element_size() * sz)
      {
        std::fill(this->begin(), this->end(), x);
      }

      template <typename OtherElementType>
      shared_base(const OtherElementType* first, const OtherElementType* last)
        : m_handle(element_size() * (last - first))
      {
        copy_typeconv(first, last, this->begin());
      }

      explicit
      shared_base(const handle_type& handle)
        : m_handle(handle)
      {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
      }

      ~shared_base() {
        if (m_handle.use_count() == 1) {
          resize(0);
        }
      }

      size_type size() const { return m_handle.size() / element_size(); }
      bool empty() const { if (size() == 0) return true; return false; }
      size_type capacity() const {
        return m_handle.capacity() / element_size();
      }
      long use_count() const { return m_handle.use_count(); }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
        reinterpret_cast<ElementType*>(m_handle.begin()), size())

            handle_type& handle()       {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        return m_handle;
      }
      const handle_type& handle() const {
        CCTBX_ARRAY_FAMILY_STATIC_ASSERT_HAS_TRIVIAL_DESTRUCTOR
        return m_handle;
      }

      shared_base<ElementType>
      deep_copy() const {
        shared_base<ElementType> result(size());
        std::copy(this->begin(), this->end(), result.begin());
        return result;
      }

      template <typename OtherElementType>
      void
      copy_from(const shared_base<OtherElementType>& a) {
        resize(a.size());
        copy_typeconv(a.begin(), a.end(), this->begin());
      }

      CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), size())

      void reserve(const size_type& new_size) {
        m_handle.reserve(element_size() * new_size);
      }

      void resize(const size_type& new_size,
                  const ElementType& x = ElementType()) {
        size_type old_size = this->size();
        if (new_size < old_size) {
          detail::destroy_array_elements(this->begin()+new_size, this->end());
        }
        m_handle.resize(element_size() * new_size);
        if (new_size > old_size) {
          std::fill(this->begin()+old_size, this->end(), x);
        }
      }

    protected:
      handle_type m_handle;
  };

}} //namespace cctbx::af

#endif //CCTBX_ARRAY_FAMILY_SHARED_BASE_H
