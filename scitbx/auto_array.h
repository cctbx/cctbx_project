#ifndef SCITBX_AUTO_ARRAY_H
#define SCITBX_AUTO_ARRAY_H

// Modified (Ralf W. Grosse-Kunstleve) copy of boost/scoped_array.hpp
// Original copyright:

//  (C) Copyright Greg Colvin and Beman Dawes 1998, 1999.
//  Copyright (c) 2001, 2002 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  http://www.boost.org/libs/smart_ptr/scoped_array.htm

#include <boost/config.hpp> // in case ptrdiff_t not in std
#include <boost/detail/workaround.hpp>
#include <cstddef> // for std::ptrdiff_t

namespace scitbx {

  //! Like std::auto_ptr, but with delete[].
  template<typename T>
  class auto_array
  {
    private:
      typedef auto_array<T> this_type;

    protected:
      mutable T* ptr;

    public:
      typedef T element_type;

      explicit
      auto_array(T* p=0) : ptr(p) {}

      // This should be non-const, but that would require serious gymnastics.
      auto_array(auto_array const& other)
      :
        ptr(const_cast<auto_array*>(&other)->release())
      {}

      // This should be non-const, but that would require serious gymnastics.
      auto_array&
      operator=(auto_array const& other)
      {
        reset(const_cast<auto_array*>(&other)->release());
        return *this;
      }

      ~auto_array() { delete[] ptr; }

      void
      reset(T* p=0) { if (p != ptr) this_type(p).swap(*this); }

      T*
      release()
      {
        T* result = ptr;
        ptr = 0;
        return result;
      }

      T&
      operator[](std::ptrdiff_t i) const { return ptr[i]; }

      T*
      get() const { return ptr; }

#if defined(__SUNPRO_CC) && BOOST_WORKAROUND(__SUNPRO_CC, <= 0x530)
      operator bool () const { return ptr != 0; }
#elif defined(__MWERKS__) \
      && BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3003))
      typedef T* (this_type::*unspecified_bool_type)() const;

      operator unspecified_bool_type() const
      {
        return ptr == 0? 0: &this_type::get;
      }
#else
      typedef T* this_type::*unspecified_bool_type;

      operator unspecified_bool_type() const
      {
        return ptr == 0? 0: &this_type::ptr;
      }

#endif

      bool
      operator!() const { return ptr == 0; }

      void swap(auto_array & other)
      {
        T* tmp = other.ptr;
        other.ptr = ptr;
        ptr = tmp;
      }
  };

  template<typename T>
  inline void
  swap(auto_array<T>& a, auto_array<T>& b) { a.swap(b); }

} // namespace scitbx

#endif  // SCITBX_AUTO_ARRAY_H
