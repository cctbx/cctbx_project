/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_VECTOR_SUB_H
#define CCTBX_VECTOR_SUB_H

#include <cassert>
#include <algorithm>

namespace cctbx {
  namespace array_family {
    namespace detail {
    /* requirements
      1. templated type
      2. explicit constructor with std::size_t
      3. allocate space & deallocate in destructor
      4. constructor from begin & end iterators
      5. begin() iterator
      6. size()
      7. resize()
      8. swap()
      9. in addition to std features, have const & non-const vecref ref()
     10. in addition to interface requirements, perform fast for this app
     11. push_back()
    */
    template <class T>
    class vector {
      public:
        typedef std::size_t size_type;
        
        explicit vector(const size_type& sz = 0) :
          m_capacity(sz+10),m_size(sz),m_data(new T[m_capacity]){}

        //constructor from begin & end iterators with copy semantics
        //don't know how this is done in std library; assume addressable 
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
    };
    } //namespace
  } //namespace
} //namespace

#endif // CCTBX_VECTOR_SUB_H
