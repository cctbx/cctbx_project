// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_H
#define CCTBX_ARRAY_H

#include <boost/array.hpp>

namespace cctbx {

  //! A cctbx::array IS-A boost::array with convenience constructors.
  template <typename T, std::size_t N>
  struct array : boost::array<T, N> {

    //! Default constructor.
    array() {}

    //! Copy boost::array with type conversion.
    template <typename U>
    array(const boost::array<U, N>& a) {
      for(std::size_t i=0;i<a.size();i++) elems[i] = a[i];
    }

    //! Convenience constructor.
    array(const T& v0
         ) {
      elems[0] = v0;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1
         ) {
      elems[0] = v0;
      elems[1] = v1;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2,
          const T& v3
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2,
          const T& v3,
          const T& v4
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
      elems[4] = v4;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2,
          const T& v3,
          const T& v4,
          const T& v5
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
      elems[4] = v4;
      elems[5] = v5;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2,
          const T& v3,
          const T& v4,
          const T& v5,
          const T& v6
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
      elems[4] = v4;
      elems[5] = v5;
      elems[6] = v6;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2,
          const T& v3,
          const T& v4,
          const T& v5,
          const T& v6,
          const T& v7
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
      elems[4] = v4;
      elems[5] = v5;
      elems[6] = v6;
      elems[7] = v7;
    }
    //! Convenience constructor.
    array(const T& v0,
          const T& v1,
          const T& v2,
          const T& v3,
          const T& v4,
          const T& v5,
          const T& v6,
          const T& v7,
          const T& v8
         ) {
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
      elems[4] = v4;
      elems[5] = v5;
      elems[6] = v6;
      elems[7] = v7;
      elems[8] = v8;
    }
    //! Convenience constructor.
    array(const T& v0,
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
      elems[0] = v0;
      elems[1] = v1;
      elems[2] = v2;
      elems[3] = v3;
      elems[4] = v4;
      elems[5] = v5;
      elems[6] = v6;
      elems[7] = v7;
      elems[8] = v8;
      elems[9] = v9;
    }
  };
}

#endif // CCTBX_ARRAY_H
