/* Copyright (c) 2002 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (N.K. Sauter)
 */

#ifndef CCTBX_SHARED_STORAGE_MEM_H
#define CCTBX_SHARED_STORAGE_MEM_H

#include <cstddef>
#include <cctbx/basic/vector.h>
#include <vector>
#include <algorithm>

namespace cctbx {

  template <typename ElementType>
  class SharedStorageMemHandle {

  public:
    typedef ElementType element_type;
    typedef std::size_t size_type;
    typedef cctbx::std_emulator::vector<element_type> vector_type;

    struct Courier {
      element_type*        px;
      size_type            sz;
      inline Courier (element_type* x, size_type s) :
        px(x),sz(s) {}
    };

    explicit SharedStorageMemHandle(const size_type& s) {
        px = new vector_type(s);
      try { 
        pn = new long(1);
      } catch (...) { delete px; throw; } 
    }  // insure no leaks if any new operator throws

    //copy constructor with reference semantics
    SharedStorageMemHandle(const SharedStorageMemHandle& r) : 
      px(r.px) { 
      ++*(pn = r.pn); 
    }

    //assignment operator with reference semantics
    SharedStorageMemHandle& operator=(const SharedStorageMemHandle& r) {
      if (pn != r.pn) {
        dispose();
        px = r.px;
        ++*(pn = r.pn);
      }
      return *this;
    }

    //constructor from a pointer and size; assume ownership of data
    //since std::vector has no constructor from a T*, I need to use
    //copy semantics, followed by destruction of the original data.
    //  ...now that I have my own emulator::vector,  I could 
    //  provide a constructor from a pointer & size
    SharedStorageMemHandle(const Courier& r) {
      //copy constructor of vector
      px = new vector_type(r.px, r.px+r.sz);
      delete r.px;
      try { 
        pn = new long(1);
      } catch (...) { delete px; throw; } 
    } 

    ~SharedStorageMemHandle() {  dispose(); }


    element_type* get() const { return &(*(px->begin())); } 
    element_type& operator[](const size_type& i) const 
                              { return (&(*(px->begin())))[i]; }

    size_type size() const            { return px->size(); }
    size_type capacity() const            { return px->capacity(); }
    long use_count() const             { return *pn; } 
    bool unique() const                { return *pn == 1; } 

    void swap(SharedStorageMemHandle<element_type>& other)  // never throws
     { std::swap(px,other.px); 
       std::swap(pn,other.pn);}

    void resize(const size_type& s)
      {px->resize(s);}

    //relinquish... and appropriate... used to transfer ownership of the       
    // referenced data from source array group to a target group of arrays 
    // all of which then share the referenced data. 
    // Note! returning a std::vector has copy semantics.  This is inefficient
    // but it's likely that I will never use the appropriate...relinquish 
    // mechanism so this can be deprecated. 

    vector_type relinquish() {
      vector_type m(1);
      px->swap(m);
      return m;
    }

    void appropriate(vector_type m) {
      px->swap(m);
    }

  protected:
    vector_type*               px;     // pointer to a std::vector
    long*                      pn;     // pointer to reference counter

    void dispose() { 
      --*pn;
      if (*pn == 0) { delete px; delete pn; }
    }


  };  //SharedStorageMemHandle

}

#endif // CCTBX_SHARED_STORAGE_MEM_H
