/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jul: created, based on a fragment from clipper/lib/hkl_datatypes.h:
               (C) 2000-2002 Kevin Cowtan
               This code is provided as free software under the CCP4
               license part (i)
 */

#ifndef CCTBX_HENDRICKSON_LATTMAN_H
#define CCTBX_HENDRICKSON_LATTMAN_H

#include <scitbx/array_family/tiny.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx {

  //! Group of Hendrickson-Lattman coefficients.
  /*! Reference: Jan Drenth,
                 Principles of Protein X-Ray Crystallography,
                 Second edition, 1999, Chapter 14.
   */
  template<typename FloatType = double>
  class hendrickson_lattman : public af::tiny_plain<FloatType, 4>
  {
    public:
      typedef af::tiny_plain<FloatType, 4> base_type;

      //! Default constructor. The coefficients are not initialized!
      hendrickson_lattman() {}

      //! Initializtion from tiny array.
      hendrickson_lattman(base_type const& coeff)
      : base_type(coeff)
      {}

      //! Initializtion from plain pointer to coefficients.
      hendrickson_lattman(const FloatType* coeff)
      {
        std::copy(coeff, coeff + 4, this->begin());
      }

      //! Construct given individual coefficients.
      hendrickson_lattman(
        FloatType const& a,
        FloatType const& b,
        FloatType const& c,
        FloatType const& d)
      {
        (*this)[0] = a;
        (*this)[1] = b;
        (*this)[2] = c;
        (*this)[3] = d;
      }

      //! Coefficients a,b,c,d as array.
      base_type const&
      coeff() const { return *this; }

      //! Coefficients a,b,c,d as array.
      base_type&
      coeff()       { return *this; }

      //! Individual coefficient a.
      FloatType const& a() const { return (*this)[0]; }

      //! Individual coefficient b.
      FloatType const& b() const { return (*this)[1]; }

      //! Individual coefficient c.
      FloatType const& c() const { return (*this)[2]; }

      //! Individual coefficient d.
      FloatType const& d() const { return (*this)[3]; }

      /*! \brief Coefficients for Friedel opposite (similar to conjugate
          complex of structure factor).
       */
      /*! Formula used: a, -b, c, -d
          <p>
          See also: cctbx::miller::sym_equiv_index
       */
      hendrickson_lattman
      conj() const
      {
        return hendrickson_lattman(a(), -b(), c(), -d());
      }

      //! Coefficients for symmetrically equivalent reflections.
      /*! The phase shift delta_phi must be given in radians.
          <p>
          See also: cctbx::miller::sym_equiv_index
       */
      hendrickson_lattman
      shift_phase(FloatType const& delta_phi) const;
  };

  template<typename FloatType>
  hendrickson_lattman<FloatType>
  hendrickson_lattman<FloatType>
  ::shift_phase(FloatType const& delta_phi) const
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

} // namespace cctbx

#endif // CCTBX_HENDRICKSON_LATTMAN_H
