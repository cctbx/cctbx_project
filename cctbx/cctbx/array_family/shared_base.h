// $Id$
/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_BASE_H
#define CCTBX_ARRAY_FAMILY_SHARED_BASE_H

#include <algorithm>
#include <cctbx/array_family/ref.h>

namespace cctbx { namespace af {
  namespace detail {

    template <typename ElementType,
              std::size_t SmallestAutoCapacity = 8>
    class basic_storage {
      public:
        typedef std::size_t size_type;

        long use_count;

        basic_storage()
          : m_capacity(0), m_size(0), m_data(0)
        {}

        explicit
        basic_storage(const size_type& sz)
          : m_capacity(sz), m_size(sz), m_data(new ElementType[sz])
        {}

        ~basic_storage() {
          delete[] m_data;
        }

              ElementType* begin()       { return m_data; }
        const ElementType* begin() const { return m_data; }

        size_type size() const { return m_size; }
        size_type capacity() const { return m_capacity; }

        void resize(const size_type& sz) {
          expand_capacity(sz);
          m_size = sz;
        }

        void auto_resize(const size_type& sz) {
          if (sz > 0 && m_capacity == 0) {
            expand_capacity(SmallestAutoCapacity);
          }
          else {
            expand_capacity(std::max(sz, m_capacity * 2));
          }
          m_size = sz;
        }

        void swap (basic_storage<ElementType>& other){
          std::swap(*this, other);
        }

      private:
        void expand_capacity(const size_type& new_capacity) {
          if (new_capacity > m_capacity) {
            m_capacity = new_capacity;
            ElementType* new_data = new ElementType[m_capacity];
            std::copy(m_data, m_data + m_size, new_data);
            delete [] m_data;
            m_data = new_data;
          }
        }

        size_type m_capacity;
        size_type m_size;
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
          : m_ptr_basic_storage(new basic_storage_type(sz)) {
          m_ptr_basic_storage->use_count = 1;
        }

        managed_storage(const managed_storage& other)
          : m_ptr_basic_storage(other.m_ptr_basic_storage) {
          m_ptr_basic_storage->use_count++;
        }

        managed_storage& operator=(const managed_storage& other) {
          if (m_ptr_basic_storage != other.m_ptr_basic_storage) {
            dispose();
            m_ptr_basic_storage = other.m_ptr_basic_storage;
            m_ptr_basic_storage->use_count++;
          }
          return *this;
        }

        ~managed_storage() { dispose(); }

        element_type* begin() const { return m_ptr_basic_storage->begin(); }

        size_type size() const { return m_ptr_basic_storage->size(); }
        size_type capacity() const { return m_ptr_basic_storage->capacity(); }
        long use_count() const { return m_ptr_basic_storage->use_count; }
        bool unique() const { return m_ptr_basic_storage->use_count == 1; }

        void swap(managed_storage<element_type>& other) {
          std::swap(m_ptr_basic_storage, other.m_ptr_basic_storage);
          std::swap(      m_ptr_basic_storage->use_count,
                    other.m_ptr_basic_storage->use_count);
        }

        void resize(const size_type& sz) {
          m_ptr_basic_storage->resize(sz);
        }

        void auto_resize(const size_type& sz) {
          m_ptr_basic_storage->resize(sz);
        }

      private:
        basic_storage_type* m_ptr_basic_storage;

        void dispose() {
          m_ptr_basic_storage->use_count--;
          if (m_ptr_basic_storage->use_count == 0) {
            delete m_ptr_basic_storage;
            m_ptr_basic_storage = 0;
          }
        }
    };

    typedef managed_storage<char> char_block;

  } //namespace detail

  template <typename ElementType>
  class shared_base
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      typedef detail::char_block handle_type;

      static size_type element_size() { return sizeof(ElementType); }

      explicit shared_base(const size_type& sz = 0)
        : m_handle( element_size() * sz )      {}

      explicit shared_base(const handle_type& handle)
        : m_handle(handle){}

      size_type size() const { return m_handle.size()/element_size(); }

      CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
        reinterpret_cast<ElementType*>(m_handle.begin()), size())

            handle_type& handle()       { return m_handle; }
      const handle_type& handle() const { return m_handle; }

      shared_base<ElementType>&
      fill(const ElementType& x) {
        std::fill(begin(), end(), x);
        return *this;
      }

      shared_base<ElementType>
      deepcopy() const {
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

      void resize(const size_type& sz) {
        m_handle.resize(element_size() * sz);
      }

    protected:
      handle_type m_handle;
  };

}} //namespace cctbx::af

#endif //CCTBX_ARRAY_FAMILY_SHARED_BASE_H
