// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_CARRAY_H
#define CCTBX_CARRAY_H

#include <boost/array.hpp>
#include <cctbx/vector/simple_io.h>

namespace cctbx {

  //! A cctbx::carray IS-A boost::array with convenience constructors.
  template <typename T, std::size_t N>
  struct carray : boost::array<T, N> {

    //! Default constructor.
    carray() {}

    //! Copy boost::array with type conversion.
    template <typename U>
    carray(const boost::array<U, N>& a) {
      for(std::size_t i=0;i<this->size();i++) this->elems[i] = a[i];
    }

    //! Convenience constructor.
    carray(const T& v0
         ) {
      this->elems[0] = v0;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3,
           const T& v4
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
      this->elems[4] = v4;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3,
           const T& v4,
           const T& v5
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
      this->elems[4] = v4;
      this->elems[5] = v5;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3,
           const T& v4,
           const T& v5,
           const T& v6
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
      this->elems[4] = v4;
      this->elems[5] = v5;
      this->elems[6] = v6;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3,
           const T& v4,
           const T& v5,
           const T& v6,
           const T& v7
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
      this->elems[4] = v4;
      this->elems[5] = v5;
      this->elems[6] = v6;
      this->elems[7] = v7;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3,
           const T& v4,
           const T& v5,
           const T& v6,
           const T& v7,
           const T& v8
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
      this->elems[4] = v4;
      this->elems[5] = v5;
      this->elems[6] = v6;
      this->elems[7] = v7;
      this->elems[8] = v8;
    }
    //! Convenience constructor.
    carray(const T& v0,
           const T& v1,
           const T& v2,
           const T& v3,
           const T& v4,
           const T& v5,
           const T& v6,
           const T& v7,
           const T& v8,
           const T& v9
         ) {
      this->elems[0] = v0;
      this->elems[1] = v1;
      this->elems[2] = v2;
      this->elems[3] = v3;
      this->elems[4] = v4;
      this->elems[5] = v5;
      this->elems[6] = v6;
      this->elems[7] = v7;
      this->elems[8] = v8;
      this->elems[9] = v9;
    }
  };

  typedef carray<int, 3> int3;
  typedef carray<int, 9> int9;
  typedef carray<double, 3> double3;
  typedef carray<double, 6> double6;
  typedef carray<double, 9> double9;

  template<class T, std::size_t N>
  std::ostream& operator<<(std::ostream& os, const carray<T, N>& v) {
    return vector::simple_output(os, v);
  }

} // namespace cctbx

#endif // CCTBX_CARRAY_H
