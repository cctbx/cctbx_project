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

namespace cctbx {
  namespace af {
    namespace detail {

    //! Partial implementation of std::vector interface.
    /*! This class is a partial implementation of std::vector with the
        extra proviso that vector<char>::begin() is correctly aligned for
        any type (same guarantees as provided by new char[]).
        Implemented interface includes:
        1. templated type
        2. explicit constructor with std::size_t
        3. allocate space in constructor & deallocate in destructor
        4. constructor from begin & end iterators
        5. begin() iterator
        6. size()
        7. resize()
        8. swap()
        9. push_back()
     */
    template <class T>
    class vector {
      public:
        typedef std::size_t size_type;

        // XXX no allocation if sz == 0
        // XXX allocate exactly sz
        explicit vector(const size_type& sz = 0) :
          m_capacity(sz+10),m_size(sz),m_data(new T[m_capacity]){}

        //!constructor from begin & end iterators with copy semantics
        // XXX allocate exactly m_size
        vector (const T* begin, const T* end) :
          m_size((end - begin)/sizeof(T)),
          m_capacity(m_size+10),
          m_data(new T[m_capacity])
          {
            std::copy(begin, end, this->begin());
          }

        ~vector() { delete[] m_data; }

        T* begin() {return m_data;}
        const T* begin() const {return m_data;}

        size_type size() const {return m_size;}
        size_type capacity() const {return m_capacity;}

        void resize(const size_type& sz) {
          if (sz>m_capacity) {
            // XXX force capacity to powers of 2, smallest something sizeable
            m_capacity = std::max(sz, m_capacity*2);
            T* bigger = new T[m_capacity];
            std::copy(m_data,m_data+m_size,bigger);
            delete [] m_data;
            m_data = bigger;
          }
          m_size = sz;
        }

        void swap (vector<T>& other){
          std::swap(*this, other);
        }

        void push_back(const T& value) {
          resize(m_size+1);
          begin()[m_size-1]=value;
        }

      private:
        size_type m_size;
        size_type m_capacity;
        T*        m_data;
    }; // class vector

    //! handle class is essentially a reference counted vector
    template <typename ElementType>
    class handle {

      public:
        typedef ElementType element_type;
        typedef std::size_t size_type;
        typedef cctbx::af::detail::vector<element_type> vector_type;

        // XXX this is like an af::ref<ElementType>
        struct courier {
          element_type*        px;
          size_type            sz;
          courier (element_type* x, size_type s) :
            px(x),sz(s) {}
        };

        explicit handle(const size_type& s) {
            px = new vector_type(s);
          try {
            pn = new long(1);
          } catch (...) { delete px; throw; }
        }  // insure no leaks if any new operator throws

        //!copy constructor with reference semantics
        handle(const handle& r) :
          px(r.px) {
          ++*(pn = r.pn);
        }

        //!assignment operator with reference semantics
        handle& operator=(const handle& r) {
          if (pn != r.pn) {
            dispose();
            px = r.px;
            ++*(pn = r.pn);
          }
          return *this;
        }

        //!constructor from a pointer and size; assume ownership of data
        /*!since std::vector has no constructor from a T*, constructor uses
           copy semantics, followed by destruction of the original data.
         */
        // XXX we cannot destroy data held by std::vector
        // XXX in what context is this useful?
        // XXX seems like handle in a handle
        handle(const courier& r) {
          px = new vector_type(r.px, r.px+r.sz);
          delete r.px;
          try {
            pn = new long(1);
          } catch (...) { delete px; throw; }
        }

        ~handle() {  dispose(); }

        element_type* get() const { return &(*(px->begin())); }
        element_type& operator[](const size_type& i) const
                                  { return (&(*(px->begin())))[i]; }

        size_type size() const            { return px->size(); }
        size_type capacity() const            { return px->capacity(); }
        long use_count() const             { return *pn; }
        bool unique() const                { return *pn == 1; }

        //! do not swap the other handle's reference counter.
        //! it is assumed that we keep the current reference
        //! count and only exchange the data.
        void swap(handle<element_type>& other)
         { std::swap(px,other.px); }

        void resize(const size_type& s)
          {px->resize(s);}

      protected:
        vector_type*               px;     // pointer to a std::vector
        long*                      pn;     // pointer to reference counter

        void dispose() {
          --*pn;
          if (*pn == 0) { delete px; delete pn; }
        }

    };  //class handle

    typedef handle<char> char_block;

    } //namespace detail


    /*!
      Hierarchy of classes:
        vector-provides a resizable array on top of a private pointer & size
        handle-has a vector and makes it shareable
        shared_base-has a handle<char>, and allows different instances
          to regard the handle as holding different data types.  Mechanism for
          cross-type handle sharing requires low-level API.
    */
    //! Share a resizable 1-dimensional data vector & regard as different types.
    template <typename ElementType>
    class shared_base
    {
      public:
        CCTBX_ARRAY_FAMILY_TYPEDEFS

        typedef detail::char_block handle_type;

        static size_type element_size() {return sizeof(ElementType);}

        explicit shared_base(const size_type& sz = 0)
          : m_handle( element_size() * sz )      {}

        //! constructor from a char block and size...reference semantics
        //! use this to share data between a group of double arrays and
        //! a group of complex<double> arrays.
        explicit shared_base(const handle_type& handle)
          : m_handle(handle){}

        //! constructor from a courier block of ElementType...assumes ownership
        explicit
        shared_base(const typename detail::handle<ElementType>::courier& c)
          : m_handle(handle_type(handle_type::courier(
                reinterpret_cast<handle_type::element_type*>(c.px),
                element_size() * c.sz
            ))){ }

        //! default copy constructor...reference semantics
        //! default assignment operator...reference semantics

        // XXX maybe these type-dependent copy semantics is too confusing?
        // XXX the user could simply use std::copy
        //! copy constructor templated on different value type...
        /*! this has deepcopy semantics...will fail unless the r
            array is initialized first XXX ???
         */
        template <class OtherElementType>
        explicit shared_base(const shared_base<OtherElementType>& r):
          m_handle (element_size() * r.size())
        {
          std::copy(r.begin(), r.end(), this->begin());
        }

        //! assignment operator templated on OtherElementType
        // ...deepcopy semantics
        template <class OtherElementType>
        shared_base<ElementType>
        operator= (const shared_base<OtherElementType>& r) {
          resize(r.size());
          std::copy(r.begin(), r.end(), this->begin());
          return *this;
        }

        size_type size() const { return m_handle.size()/element_size(); }

        CCTBX_ARRAY_FAMILY_BEGIN_END_ETC(
          reinterpret_cast<ElementType*>(m_handle.get()), size())

              handle_type& handle()       { return m_handle; }
        const handle_type& handle() const { return m_handle; }

        shared_base<ElementType>
        deepcopy() const {
          shared_base<ElementType> result(size());
          std::copy(this->begin(), this->end(), result.begin());
          return result;
        }

        CCTBX_ARRAY_FAMILY_TAKE_REF(begin(), size())

        //! Resize all shared copies; single-instance deepcopy semantics;
        //! extra elements not initialized
        void resize(const size_type& sz = 0) {
          m_handle.resize(element_size() * sz);
        }
      protected:
        handle_type m_handle;
    };

  } //namespace af
} //namespace cctbx

#endif //CCTBX_ARRAY_FAMILY_SHARED_BASE_H
