// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: created, based on a fragment from clipper/lib/hkl_datatypes.h:
               (C) 2000-2002 Kevin Cowtan
               This code is provided as free software under the CCP4
               license part (i)
 */

#ifndef CCTBX_HENDRICKSON_LATTMAN_H
#define CCTBX_HENDRICKSON_LATTMAN_H

#include <cctbx/fixes/cmath>
#include <cctbx/array_family/tiny.h>

namespace cctbx {

  //! Group of Hendrickson-Lattman coefficients.
  /*! Reference: Jan Drenth,
                 Principles of Protein X-Ray Crystallography,
                 Second edition, 1999, Chapter 14.
   */
  template<typename FloatType>
  class hendrickson_lattman
  {
    public:
      //! Default constructor. The coefficients are not initialized!
      hendrickson_lattman() {}

      //! Construct from tiny array.
      hendrickson_lattman(af::tiny<FloatType, 4> const& coeff)
      : coeff_(coeff)
      {}

      //! Construct from plain pointer to coefficients.
      hendrickson_lattman(const FloatType* coeff)
      {
        std::copy(coeff, coeff + 4, coeff_);
      }

      //! Construct given individual coefficients.
      hendrickson_lattman(
        FloatType const& a,
        FloatType const& b,
        FloatType const& c,
        FloatType const& d)
      {
        coeff_[0] = a;
        coeff_[1] = b;
        coeff_[2] = c;
        coeff_[3] = d;
      }

      //! Coefficients a,b,c,d as array.
      af::tiny<FloatType, 4> const& array() const { return coeff_; }

      //! Individual coefficient a.
      FloatType const& a() const { return coeff_[0]; }
      //! Individual coefficient b.
      FloatType const& b() const { return coeff_[1]; }
      //! Individual coefficient c.
      FloatType const& c() const { return coeff_[2]; }
      //! Individual coefficient d.
      FloatType const& d() const { return coeff_[3]; }

      /*! \brief Coefficients for Friedel opposite (similar to conjugate
          complex of structure factor).
       */
      /*! Formula used: a, -b, c, -d
          <p>
          See also: class cctbx::Miller::SymEquivIndex
       */
      hendrickson_lattman
      conj() const
      {
        return hendrickson_lattman(a(), -b(), c(), -d());
      }

      //! Coefficients for symmetrically equivalent reflections.
      /*! The phase shift delta_phi must be given in radians.
          <p>
          See also: class cctbx::Miller::SymEquivIndex
       */
      hendrickson_lattman
      shift_phase(FloatType const& delta_phi) const
      {
        FloatType c1 = std::cos(delta_phi);
        FloatType s1 = std::sin(delta_phi);
        FloatType c2 = std::cos(2. * delta_phi);
        FloatType s2 = std::sin(2. * delta_phi);
        return hendrickson_lattman(
          a()*c1 - b()*s1,
          a()*s1 + b()*c1,
          c()*c2 - d()*s2,
          c()*s2 + d()*c2);
      }

    private:
      af::tiny<FloatType, 4> coeff_;
  };

} // namespace cctbx

#endif // CCTBX_HENDRICKSON_LATTMAN_H
