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

#include <cctbx/array_family/tiny_types.h>
#include <cctbx/array_family/tiny_algebra.h>

namespace cctbx {

  //! Class for cartesian (orthogonal, real) coordinates.
  /*! The template parameter FloatType should be a floating point type
      (e.g. float or double).
      <p>
      See also: class fractional
   */
  template <class FloatType>
  class cartesian : public af::tiny<FloatType, 3> {
    public:
      //! The elements of the coordinate vector are initialized with 0.
      cartesian() {
        for(std::size_t i=0;i<3;i++) this->elems[i] = 0;
      }
      //! The elements of the coordinate vector are copied from v.
      template <class OtherFloatType>
      cartesian(const af::tiny<OtherFloatType, 3> v) {
        for(std::size_t i=0;i<3;i++) this->elems[i] = v[i];
      }
      //! The elements of the coordinate vector are copied from xyz.
      template <class OtherFloatType>
      cartesian(const OtherFloatType* xyz) {
        for(std::size_t i=0;i<3;i++) this->elems[i] = xyz[i];
      }
      //! The elements of the coordinate vector are initialized with x,y,z.
      cartesian(const FloatType& x, const FloatType& y, const FloatType& z) {
        this->elems[0] = x; this->elems[1] = y; this->elems[2] = z;
      }
      //! Length squared (scalar product) of the coordinate vector.
      FloatType Length2() const {
        return af::sum((*this) * (*this));
      }
  };

  namespace detail {

    template <class FloatType>
    af::int3
    getUnitShifts(const af::tiny<FloatType, 3>& Delta)
    {
      af::int3 result;
      for(std::size_t i=0;i<3;i++) {
        if (Delta[i] >= 0.) result[i] = int(Delta[i] + 0.5);
        else                result[i] = int(Delta[i] - 0.5);
      }
      return result;
    }

  } // namespace detail

  //! Class for fractional coordinates.
  /*! The template parameter FloatType should be a floating point type
      (e.g. float or double).
      <p>
      See also: class cartesian
   */
  template <class FloatType>
  class fractional : public af::tiny<FloatType, 3> {
    public:
      //! The elements of the coordinate vector are initialized with 0.
      fractional() {
        for(std::size_t i=0;i<3;i++) this->elems[i] = 0;
      }
      //! The elements of the coordinate vector are copied from v.
      template <class OtherFloatType>
      fractional(const af::tiny<OtherFloatType, 3> v) {
        for(std::size_t i=0;i<3;i++) this->elems[i] = v[i];
      }
      //! The elements of the coordinate vector are copied from xyz.
      template <class OtherFloatType>
      fractional(const OtherFloatType* xyz) {
        for(std::size_t i=0;i<3;i++) this->elems[i] = xyz[i];
      }
      //! The elements of the coordinate vector are initialized with x,y,z.
      fractional(const FloatType& x, const FloatType& y, const FloatType& z) {
        this->elems[0] = x; this->elems[1] = y; this->elems[2] = z;
      }
      /*! \brief Apply modulus operation such that 0.0 <= x < 1.0
          for all elements of the coordinate vector.
       */
      fractional modPositive() const {
        fractional result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = std::fmod(this->elems[i], 1.);
          while (result[i] <  0.) result[i] += 1.;
          while (result[i] >= 1.) result[i] -= 1.;
        }
        return result;
      }
      /*! \brief Apply modulus operation such that -0.5 < x <= 0.5
          for all elements of the coordinate vector.
       */
      fractional modShort() const {
        fractional result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = std::fmod(this->elems[i], 1.);
          if      (result[i] <= -.5) result[i] += 1.;
          else if (result[i] >   .5) result[i] -= 1.;
        }
        return result;
      }
  };

} // namespace cctbx

#endif // CCTBX_COORDINATES_H
