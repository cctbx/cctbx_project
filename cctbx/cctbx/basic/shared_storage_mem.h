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
#include <vector>
#include <algorithm>

namespace cctbx {

  template <typename ElementType>
  class SharedStorageMemHandle {

  public:
    typedef ElementType element_type;
    typedef std::size_t size_type;

    struct Courier {
      element_type*        px;
      size_type            sz;
      inline Courier (element_type* x, size_type s) :
        px(x),sz(s) {}
    };

    explicit SharedStorageMemHandle(const size_type& s) {
      element_type* pi;
      //try{
        ps = new size_type(s);
      //} catch (...) { throw; }
      try{
        pi = new element_type[s];
      } catch (...) { delete ps; throw; }
      try {
        px = new element_type*(pi);
      } catch (...) { delete ps; delete [] pi; throw; }
      try { 
        pn = new long(1);
      } catch (...) { delete ps; delete [] pi; delete px; throw; } 
    }  // pedanticly insures no leaks if any new operator throws

    //copy constructor with reference semantics
    SharedStorageMemHandle(const SharedStorageMemHandle& r) : 
      px(r.px), ps(r.ps) { 
      ++*(pn = r.pn); 
    }

    //assignment operator with reference semantics
    SharedStorageMemHandle& operator=(const SharedStorageMemHandle& r) {
      if (pn != r.pn) {
        dispose();
        px = r.px;
        ps = r.ps;
        ++*(pn = r.pn);
      }
      return *this;
    }

    //constructor from a pointer and size; assume ownership of data
    SharedStorageMemHandle(const Courier& r) {
      ps = new size_type(r.sz);
      try {
        px = new element_type*(r.px);
      } catch (...) { delete ps; delete [] r.px; throw; }
      try { 
        pn = new long(1);
      } catch (...) { delete ps; delete [] *px; delete px; throw; } 
    } 
      

    ~SharedStorageMemHandle() {  dispose(); }


    element_type* get() const                 { return *px; } 
    element_type& operator[](const size_type& i) const { return (*px)[i]; }

    size_type size() const             { return *ps; }
    long use_count() const             { return *pn; } 
    bool unique() const                { return *pn == 1; } 

    void swap(SharedStorageMemHandle<element_type>& other)  // never throws
     { std::swap(px,other.px); 
       std::swap(pn,other.pn); std::swap(ps,other.ps);}

    //relinquish... and appropriate... used to transfer ownership of the       
    // referenced data from source array group to a target group of arrays 
    // all of which then share the referenced data.  

    Courier relinquish() {
      Courier m(*px,*ps);
      delete px;
      element_type* pi;
      try {
        pi = new element_type[1];
      } catch (...) { /* no real way to recover from this? */; throw; }
      *ps = 1;
      try {
        px = new element_type*(pi);
      } catch (...) { delete [] pi; throw; }
      return m;
    } // Courier has no memory management; relinquished data must be 
      // appropriated by exactly one SharedStorageMemHandle before Courier
      // goes out of scope.

    void appropriate(const Courier& m) {
      delete [] *px;
      *px = m.px;
      *ps = m.sz;
    }

  protected:

    element_type**    px;     // pointer to a pointer to elements
    long*             pn;     // pointer to reference counter
    size_type*        ps;     // pointer to the size lookup 

    void dispose() { 
      --*pn;
      if (*pn == 0) { delete [] *px; delete px; delete pn; delete ps; }
    }


  };  //SharedStorageMemHandle

}

#endif // CCTBX_SHARED_STORAGE_MEM_H
