/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_AF_DYNAMIC_H
#define CCTBX_AF_DYNAMIC_H

#include <cassert>
#include <algorithm>
#include <cstddef>

#include <cctbx/carray.h>
#include <cctbx/vector/reductions.h>
#include <cctbx/vecref.h>

namespace cctbx {
  namespace array_family {
    namespace detail {
    //! Lightweight implementation of std::vector for improved performance
    /*! This class is a partial implementation of std::vector with the 
        extra proviso that vector<char> performs at the same speed as 
        char* for the application.  Implemented interface includes:
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
        
        explicit vector(const size_type& sz = 0) :
          m_capacity(sz+10),m_size(sz),m_data(new T[m_capacity]){}

        //!constructor from begin & end iterators with copy semantics
        vector (const T* begin, const T* end) :
          m_size(((long)end - (long)begin)/sizeof(T)), 
          m_capacity(m_size+10),
          m_data(new T[m_capacity])
          {
            for (int i=0; i<m_size; ++i) {m_data[i]=begin[i];} 
          }

        ~vector() { delete[] m_data; }

        T* begin() {return m_data;}
        const T* begin() const {return m_data;}

        const size_type size() const {return m_size;}
        const size_type capacity() const {return m_capacity;}

        void resize(const size_type& sz) {
          assert (sz>m_size);
          if (sz>m_capacity) {
            m_capacity = std::max(sz, m_capacity*2);
            T* bigger = new T[m_capacity];
            std::copy(m_data,m_data+m_size,bigger);
            T* old = m_data;
            m_data = bigger;
            delete [] old;
          }
          m_size = sz;
        }

        void swap (vector<T>& other){
          std::swap(m_data, other.m_data);
          std::swap(m_size, other.m_size);
          std::swap(m_capacity, other.m_capacity);          
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
        typedef cctbx::array_family::detail::vector<element_type> vector_type;

        struct Courier {
          element_type*        px;
          size_type            sz;
          inline Courier (element_type* x, size_type s) :
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
        handle(const Courier& r) {
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

    } //namespace detail

    typedef detail::handle<char> char_block;

    /*!
      Hierarchy of classes:
        vector-provides a resizable array on top of a private pointer & size
        handle-has a vector and makes it shareable
        dynamic-has a handle<char>, and allows different instances
          to regard the handle as holding different data types.  Mechanism for
          cross-type handle sharing requires low-level API.  
    */
    //! Share a resizable 1-dimensional data vector & regard as different types.
    template <typename ValueType>
    class dynamic             
    {
      public:
        typedef ValueType        value_type;
        typedef ValueType*       iterator;
        typedef const ValueType* const_iterator;
        typedef ValueType&       reference;
        typedef const ValueType& const_reference;
        typedef std::size_t      sz_t;
        typedef char_block       handle_type;

        const sz_t element_size() const {return sizeof(value_type);}

        explicit dynamic(const sz_t& sz = 0)
          : m_handle( element_size() * sz )      {      }

        //! constructor from a char block and size...reference semantics
        //! use this to share data between a group of double arrays and 
        //! a group of complex<double> arrays.  
        dynamic(const handle_type& handle)
          : m_handle(handle){}

        //! constructor from a Courier block of value_type...assumes ownership
        dynamic(const typename detail::handle<ValueType>::Courier& c)
          : m_handle(handle_type(handle_type::Courier(
                reinterpret_cast<handle_type::element_type*>(c.px),
                element_size() * c.sz
            ))){ }

        //! default copy constructor...reference semantics 
        //! default assignment operator...reference semantics

        //! copy constructor templated on different value type...
        /*! this has deepcopy semantics...will fail unless the r
            array is initialized first
         */
        template <class ValueType2>
        explicit dynamic(const dynamic<ValueType2>& r):
          m_handle (element_size() * r.size())
        {
          std::cout<<"...templated copy constructor"<<std::endl;
          std::copy(r.begin(), r.end(), this->begin());
        }

        //! assignment operator templated on value_type2...deepcopy semantics
        //! insist that lhs be fully constructed first with appropriate size()
        //! in general rhs must be initialized first
        template <class ValueType2>
        dynamic<value_type> operator= (const dynamic<ValueType2>& r) {
          std::cout<<"...templated assignment operator"<<std::endl;
          assert (this->size() == r.size());        
          std::copy(r.begin(), r.end(), this->begin());
          return *this;
        }

              value_type* begin()     
                { return reinterpret_cast<value_type*>(m_handle.get());}
        const value_type* begin() const 
                { return reinterpret_cast<value_type*>(m_handle.get()); }
              value_type* end()       { return begin() + size(); }
        const value_type* end() const { return begin() + size(); }

              value_type& operator[](sz_t i)       { return begin()[i]; }
        const value_type& operator[](sz_t i) const { return begin()[i]; }

        sz_t size() const { return m_handle.size()/element_size(); }

              handle_type& handle()       { return m_handle; }
        const handle_type& handle() const { return m_handle; }

        //! Use vecref for faster performance; mimics raw pointers
        vecref<value_type> ref()       {
          return vecref<      value_type>(begin(), size());
        }
        vecref<const value_type> ref() const {
          return vecref<const value_type>(begin(), size());
        }

        dynamic<value_type>
        deepcopy() const
        {
          dynamic<value_type> result(size());
          std::copy(this->begin(), this->end(), result.begin());
          return result;
        }

        //! Resize all shared copies; single-instance deepcopy semantics; 
        //! extra elements not initialized
        void resize(const sz_t& sz = 0) {
          assert (sz > size()); // permit expansion only
          m_handle.resize(element_size() * sz); 
        }

        void push_back(const value_type& value) {
          m_handle.resize(size()+1);
          this->operator[](size()-1)=value;
        }

        //! swaps in the data from another handle but preserves reference count
        void swap(dynamic<value_type>& other) {
          m_handle.swap(other.m_handle);
        }

        //! assign a scalar value to all elements of the array
        //! default value is taken from default constructor of value_type
        void assign(const value_type& value = value_type()) {
          for (iterator position = begin(); position!=end(); ++position) {
            *position = value;
          }
        }
      protected:
        handle_type m_handle;
    };

  } //namespace array_family
} //namespace


#endif //CCTBX_AF_DYNAMIC_H
