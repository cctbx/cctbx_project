// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Created 2001 Jul 03 (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_COORDINATES_H
#define CCTBX_COORDINATES_H

#include <boost/array.hpp>
#include <cctbx/basic/matrixlite.h>

namespace cctbx {

  template <class T>
  class cartesian : public boost::array<T, 3> {
    public:
      cartesian() {
        for(std::size_t i=0;i<3;i++) elems[i] = 0;
      }
      template <class U>
      cartesian(const boost::array<U, 3> v) {
        for(std::size_t i=0;i<3;i++) elems[i] = v[i];
      }
      template <class U>
      cartesian(const U* xyz) {
        for(std::size_t i=0;i<3;i++) elems[i] = xyz[i];
      }
      cartesian(const T& x, const T& y, const T& z) {
        elems[0] = x; elems[1] = y; elems[2] = z;
      }
      inline T Length2() const {
        return (*this) * (*this);
      }
  };

  template <class T>
  class fractional : public boost::array<T, 3> {
    public:
      fractional() {
        for(std::size_t i=0;i<3;i++) elems[i] = 0;
      }
      template <class U>
      fractional(const boost::array<U, 3> v) {
        for(std::size_t i=0;i<3;i++) elems[i] = v[i];
      }
      template <class U>
      fractional(const U* xyz) {
        for(std::size_t i=0;i<3;i++) elems[i] = xyz[i];
      }
      fractional(const T& x, const T& y, const T& z) {
        elems[0] = x; elems[1] = y; elems[2] = z;
      }
      fractional modPositive() const {
        fractional result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = std::fmod(elems[i], 1.);
          while (result[i] <  0.) result[i] += 1.;
          while (result[i] >= 1.) result[i] -= 1.;
        }
        return result;
      }
      fractional modShort() const {
        fractional result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = std::fmod(elems[i], 1.);
          if      (result[i] <= -.5) result[i] += 1.;
          else if (result[i] >   .5) result[i] -= 1.;
        }
        return result;
      }
  };

} // namespace cctbx

#endif // CCTBX_COORDINATES_H
