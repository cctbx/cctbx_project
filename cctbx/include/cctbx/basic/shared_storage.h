/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_SHARED_STORAGE_H
#define CCTBX_SHARED_STORAGE_H

#include <cstddef>
#include <vector>
#include <algorithm>
#include <cassert> //for now--Later use cctbx_assert instead
#include <iostream> //for now

#include <cctbx/basic/shared_storage_mem.h>
#include <cctbx/carray.h>
#include <cctbx/vector/reductions.h>
#include <cctbx/vecref.h>

namespace cctbx {

  typedef SharedStorageMemHandle<char> char_block;

  /*
    Hierarchy of classes:
      vector-provides a resizable array on top of a private pointer & size
      SharedStorageMemHandle-has a vector, and makes it shareable
      SharedStorage-has a SharedStorageMemHandle, and allows different instances
        to regard the handle as holding different data types.  Mechanism for
        cross-type handle sharing requires low-level API.  
      SharedStorageND-is a SharedStorage, but also has a shared N-dimensional
        shape (extent) parameter defining how the N-dimensional array is 
        layed out in memory.
  */

  template <typename ValueType>
  class SharedStorage             // A one-dimensional shared array
  {
    public:
      typedef ValueType        value_type;
      typedef ValueType*       iterator;
      typedef const ValueType* const_iterator;
      typedef ValueType&       reference;
      typedef const ValueType& const_reference;
      typedef std::size_t      sz_t;
      typedef std::ptrdiff_t   difference_type;
      typedef char_block       handle_type;

      const sz_t element_size() const {return sizeof(value_type);}

      explicit SharedStorage(const sz_t& sz = 0)
        : m_handle( element_size() * sz )      {      }

      // constructor from a char block and size...reference semantics
      // use this to share data between a group of double arrays and a group of
      // complex<double> arrays.  Note: The use of sz for defining a private size
      // has been deprecated.  Slicing now requires a special View-class 
      SharedStorage(const handle_type& handle)
        : m_handle(handle){}

      // constructor from a Courier block of value_type...assumes ownership
      SharedStorage(const typename SharedStorageMemHandle<ValueType>::Courier& c)
        : m_handle(handle_type(handle_type::Courier(
              reinterpret_cast<handle_type::element_type*>(c.px),
              element_size() * c.sz
          ))){ }

      // default copy constructor...reference semantics 
      // default assignment operator...reference semantics

      // copy constructor templated on different value type...deepcopy semantics
      // don't try this without initializing r first
      template <class ValueType2>
      explicit SharedStorage(const SharedStorage<ValueType2>& r):
        m_handle (element_size() * r.size())
      {
        std::copy(r.begin(), r.end(), this->begin());
      }

      // assignment operator templated on value_type2...deepcopy semantics
      // insist that lhs be fully constructed first with appropriate size()
      // in general rhs must be initialized first
      template <class ValueType2>
      SharedStorage<value_type> operator= (const SharedStorage<ValueType2>& r) {
        std::cout<<"in templated op="<<std::endl;
        assert (this->size() == r.size());        
        std::copy(r.begin(), r.end(), this->begin());
        return *this;
      }


      // Put this back in later
      //template <typename T>
      //shared_storage(const std::vector<T>& std_vec)
      //  : m_handle(new char[element_size() * std_vec.size()]),
      //    m_size(std_vec.size())
      //{
      //  m_begin = reinterpret_cast<value_type*>(m_handle.get());
      //  std::copy(std_vec.begin(), std_vec.end(), m_begin);
      //}

            value_type* begin()       //{return m_begin;}
              { return reinterpret_cast<value_type*>(m_handle.get());}
      const value_type* begin() const //{return m_begin;}
              { return reinterpret_cast<value_type*>(m_handle.get()); }
            value_type* end()       { return begin() + size(); }
      const value_type* end() const { return begin() + size(); }

            value_type& operator[](sz_t i)       { return begin()[i]; }
      const value_type& operator[](sz_t i) const { return begin()[i]; }
      

      //semantics:  SharedStorage is a 1-D shared array in which the
      //size is completely determined by the underlying mem_handle.  If a subset
      //is desired, a new View class must be instantiated that has a private size.  
      sz_t size() const { return m_handle.size()/element_size(); }

            handle_type& handle()       { return m_handle; }
      const handle_type& handle() const { return m_handle; }

      // Use vecref for faster performance; mimics raw pointers
      vecref<value_type> ref()       {
        return vecref<      value_type>(begin(), size());
      }
      vecref<const value_type> ref() const {
        return vecref<const value_type>(begin(), size());
      }

      SharedStorage<value_type>
      deepcopy() const
      {
        SharedStorage<value_type> result(size());
        std::copy(this->begin(), this->end(), result.begin());
        return result;
      }

      // Resize all shared copies; single-instance deepcopy semantics; 
      // extra elements not initialized
      void resize(const sz_t& sz = 0) {
        assert (sz > size()); // permit expansion only
        m_handle.resize(element_size() * sz); 
      }

      void push_back(const value_type& value) {
        m_handle.resize(size()+1);
        this->operator[](size()-1)=value;
      }

    protected:
      handle_type m_handle;
  };


  template <typename ValueType, int D>
  class SharedStorageND: public SharedStorage<ValueType>
  {
    public:
      typedef SharedStorageND<ValueType,D>        SSND_type;
      typedef typename SSND_type::SharedStorage   SS_type;
      typedef typename 
       SharedStorageMemHandle<ValueType>::Courier courier_type;
      typedef std::size_t                         sz_t;
      typedef SharedStorageMemHandle<int>         shared_tuple_type;
      typedef carray<int, D>                      external_tuple_type;
      typedef ValueType                           value_type;
      typedef ValueType*                          iterator;
      typedef const ValueType*                    const_iterator;
      typedef ValueType&                          reference;
      typedef const ValueType&                    const_reference;
      typedef char_block                          handle_type;
      typedef SharedStorageMemHandle<int>         sz_handle_type;

      // element_size() - inherited

      SharedStorageND()
        : SS_type( 0 ),
          m_shape(D)
      {           for (int i = 0; i<D; ++i) {m_shape[i]=0;}   }

      explicit SharedStorageND(const sz_t& sz)
        : SS_type( sz ),
          m_shape(D)
      {  assert (D==1); m_shape[0]=sz * this->element_size();   }

      explicit SharedStorageND(const sz_t& sz0, const sz_t& sz1)
        : SS_type( sz0 * sz1 ),
          m_shape(D)
      {  assert (D==2); 
         m_shape[0]=sz0;  
         m_shape[1]=sz1 * this->element_size();  
      }

      explicit SharedStorageND(const sz_t& sz0, const sz_t& sz1, const sz_t& sz2)
        : SS_type( sz0 * sz1 * sz2 ),
          m_shape(3)
      {  assert (D==3);  
         m_shape[0]=sz0;  
         m_shape[1]=sz1;  
         m_shape[2]=sz2 * this->element_size();  
      }

      SharedStorageND(const external_tuple_type& t)
        : SS_type( cctbx::vector::product(t) ),
          m_shape(D)
      {
         for (int i = 0; i<D; ++i) {m_shape[i]=t[i];}
         m_shape[D-1]*=this->element_size();
      }

      // constructor from a shared char block and shared N-dim size
      // reference semantics
      // This constructor is the reason for having m_shape as private 
      // data...it provides storage extent information shared even among instances
      // having different value_type 
      SharedStorageND(const handle_type& handle, const shared_tuple_type& sz)
        : SS_type(handle),
          m_shape(sz)
      {
        assert (cctbx::vector::product(this->shape()) * this->element_size() 
                 == this->m_handle.size());
      }

      // constructor from a Courier block of value_type...assumes ownership
      SharedStorageND(const courier_type& c, 
                      const external_tuple_type& sz)
        : SS_type(handle_type(handle_type::Courier(
              reinterpret_cast<handle_type::element_type*>(c.px),
              this->element_size() * cctbx::vector::product(sz)
          ))), m_shape(D){
          for (int i = 0; i<D; ++i) m_shape[i]=sz[i];
          m_shape[D-1]*=this->element_size();
      }
 
      // default copy constructor...reference semantics 
      // default assignment operator...reference semantics

      // copy constructor templated on different value type...deepcopy semantics
      // don't try this without initializing r first
      template <class ValueType2>
      explicit SharedStorageND(const SharedStorageND<ValueType2,D>& r):
        SS_type (cctbx::vector::product(r.shape())),
        m_shape(D)
      {
        for (int i = 0; i<D; ++i) m_shape[i]=r.shape()[i];
        m_shape[D-1]*=this->element_size();
        std::copy(r.begin(), r.end(), this->begin());
      }

      // assignment operator templated on value_type2...deepcopy semantics
      // insist that lhs be fully constructed first with appropriate size()
      // in general rhs must be initialized first
      template <class ValueType2>
      SSND_type operator=(const SharedStorageND<ValueType2,D>& r) 
      {
        std::cout<<"in templated op="<<std::endl;
        m_shape = shared_tuple_type(D);
        for (int i = 0; i<D; ++i) {m_shape[i]=r.shape()[i];}
        m_shape[D-1]*=this->element_size();
        this->m_handle= handle_type( this->element_size() * 
                               cctbx::vector::product(this->shape()));
        std::copy(r.begin(), r.end(), this->begin());
        return *this;
      }

      sz_t size() const { 
            return cctbx::vector::product(this->shape()); }

      // shape() returns an external_tuple_type that is effectively
      // the integer dimensions of the array; e.g. at 3x4x6 array
      const external_tuple_type shape() const {
        external_tuple_type result;
        for (int i =0; i<D; ++i) {result[i] = m_shape[i];}
        result[D-1]/=this->element_size();
        return result;
      }

      // size_handle() returns the shared dimensions of the underlying
      // character block.  This is necessary when sharing data between 
      // two arrays of different type.  
      // For example, suppose a (double) array has shape() = 3x4x6
      // A (complex<double>) shared array would have shape() = 3x4x3
      // The size_handle() for both arrays would be 3x4x48
         sz_handle_type& size_handle()       { return m_shape; }
   const sz_handle_type& size_handle() const { return m_shape; }


      SSND_type
      deepcopy() const
      {
        external_tuple_type extent = this->shape();
        SSND_type result(extent);
        std::copy(this->begin(), this->end(), result.begin());
        return result;
      }

      // Resize all shared copies; single-instance deepcopy semantics; 
      // extra elements initialized!!!
      // This function is the whole raison d'etre of the special
      // SharedStorageND class.  
/*
      void resize(const sz_t& sz = 0) {
        assert (sz > size()); // permit expansion only
        SharedStorage<value_type> newstorage(sz);
        std::copy(this->begin(), this->end(), newstorage.begin());
        m_handle.appropriate(newstorage.handle().relinquish());
        m_size = sz;
      }
*/

    protected:
      shared_tuple_type  m_shape;
  };

 

}

#endif // CCTBX_SHARED_STORAGE_H
