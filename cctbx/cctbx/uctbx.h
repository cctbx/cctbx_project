// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for the handling of unit cell parameters.
 */

#ifndef CCTBX_UCTBX_H
#define CCTBX_UCTBX_H

#include <iostream>
#include <cctbx/fixes/cmath>
#include <boost/array.hpp>
#include <cctbx/error.h>
#include <cctbx/basic/matrixlite.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/matrix.h>
#include <cctbx/constants.h>

//! Unit Cell Toolbox namespace.
namespace uctbx {

  using namespace cctbx;

  static const error
    corrupt_unit_cell_parameters("Corrupt unit cell parameters.");
  static const error
    corrupt_metrical_matrix("Corrupt metrical matrix.");

  using MatrixLite::dtype::Vec3;
  using MatrixLite::dtype::Mx33;

  //! @name Conversions between radians and degrees.
  //@{
  inline double deg_as_rad(double deg) { return deg * constants::pi_180;}
  inline double rad_as_deg(double rad) { return rad / constants::pi_180;}
  //@}

  //! inline function for fast matrix * vector computation.
  inline Vec3 operator*(const Mx33& m, const Vec3& v)
  {
    Vec3 mv;
    mv[0] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
    mv[1] = m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
    mv[2] = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
    return mv;
  }

  //! Helper function for transforming metrical matrices.
  inline Mx33 getRtGR(const Mx33& G, const Mx33& R)
  {
    Mx33 Rt, GR, RtGR;
    MatrixLite::transpose<double>(R.elems, 3, 3, Rt.elems);
    MatrixLite::multiply<double>(G.elems, R.elems, 3, 3, 3, GR.elems);
    MatrixLite::multiply<double>(Rt.elems, GR.elems, 3, 3, 3, RtGR.elems);
    return RtGR;
  }

  //! @name inline functions for fast conversions of d-spacing measures.
  //@{
  inline double Q_as_s(double Q) {  return std::sqrt(Q); }
  inline double Q_as_d(double Q) {
    if (Q == 0.) return -1.;
    /* else */   return 1. / std::sqrt(Q);
  }
  //@}

  //! Helper class for passing unit cell parameters.
  class uc_params : public boost::array<double, 6> {
    public:
      //! @name Constructors.
      //@{
      //! Constructor using parameters (a, b, c, alpha, beta, gamma).
      uc_params(double a = 1.,
                double b = 1.,
                double c = 1.,
                double alpha = 90.,
                double beta  = 90.,
                double gamma = 90.) {
        elems[0] = a; elems[1] = b; elems[2] = c;
        elems[3] = alpha; elems[4] = beta; elems[5] = gamma;
      }
      //! Constructor using arrays for lengths and angles.
      uc_params(const Vec3& Len, const Vec3& Ang) {
        int i;
        for(i=0;i<3;i++) elems[i]     = Len[i];
        for(i=0;i<3;i++) elems[i + 3] = Ang[i];
      }
      //@}
      //! @name Access to arrays of lengths and angles.
      //@{
      //!
      inline double* Len() { return &elems[0]; }
      inline double* Ang() { return &elems[3]; }
      inline const double* Len() const { return &elems[0]; }
      inline const double* Ang() const { return &elems[3]; }
      //@}
      //! @name Access to individual lengths and angles.
      //@{
      inline double Len(int i) const { return elems[i]; }
      inline double Ang(int i) const { return elems[3 + i]; }
      //@}
  };

  //! Main class for the handling of unit cell information.
  /*! All angels accessed through the public interface are in degrees.
      Internally, the angles are stored in radians.
      <p>
      The PDB convention for orthogonalization and fractionalization
      of coordinates is used:
      <pre>
      Crystallographic Basis: D = {a,b,c}
      Cartesian Basis:        C = {i,j,k}
      i || a
      j is in (a,b) plane
      k = i x j
      </pre>
    */
  class UnitCell {

    public:
      //! @name Constructors.
      //@{
      //! Constructor using parameters (a, b, c, alpha, beta, gamma).
      UnitCell(const uc_params& ucp);
      //! Constructor using parameters derived from a metrical matrix.
      /*! The metrical matrix is defined as:
         <pre>
         ( a*a,            a*b*cos(gamma), a*c*cos(beta)  )
         ( a*b*cos(gamma), b*b,            b*c*cos(alpha) )
         ( a*c*cos(beta),  b*c*cos(alpha), c*c            )
         </pre>
       */
      UnitCell(const Mx33& MetricalMatrix);
      //@}

      //! @name Query parameters and volume.
      //@{
      uc_params getParameters(bool reciprocal = false) const;
      inline const Mx33& getMetricalMatrix(bool reciprocal = false) const {
        if (reciprocal == false) return G;
        /* else */               return R_G;
      }
      inline double getVolume() const { return Vol; }
      //@}

      //! @name Orthogonalization and fractionalization of coordinates.
      //@{
      //! This matrix converts cartesian to fractional coordinates.<br>
      //! x(fractional) = matrix * x(cartesian).
      inline const Mx33& getFractionalizationMatrix() const { return Frac; }
      //! This matrix converts fractional to cartesian coordinates.<br>
      //! x(cartesian) = matrix * x(fractional).
      inline const Mx33& getOrthogonalizationMatrix() const { return Orth; }
      //! Converts cartesian coordinates Xc to fractional coordinates.
      inline Vec3 fractionalize(const Vec3& Xc) const { return Frac * Xc; }
      //! Converts fractional coordinates Xf to cartesian coordinates.
      inline Vec3 orthogonalize(const Vec3& Xf) const { return Orth * Xf; }
      //@}

      //! @name Transformation (change-of-basis) of unit cell parameters.
      //@{
      //! InvCBMxR is the inverse of the 3x3 change-of-basis matrix
      //! that transforms coordinates in the old basis system to
      //! coodinates in the new basis system.
      UnitCell ChangeBasis(const Mx33& InvCBMxR, double RBF = 1.) const;
      UnitCell ChangeBasis(const sgtbx::RotMx& InvCBMxR) const;
      //@}

      //! @name Methods using Miller indices.
      //@{
      //! Compute the maximum Miller indices for a given minimum d-spacing.
      Miller::Index MaxMillerIndices(double dmin) const;
      //! d-spacing measure Q = 1 / d^2 = s^2 = (2*sin(theta)/lambda)^2.
      inline double Q(const Miller::Index& MIx) const
      {
        return
            (MIx[0] * MIx[0]) * (R_Len[0] * R_Len[0])
          + (MIx[1] * MIx[1]) * (R_Len[1] * R_Len[1])
          + (MIx[2] * MIx[2]) * (R_Len[2] * R_Len[2])
          + (2 * MIx[0] * MIx[1]) * (R_Len[0] * R_Len[1] * R_cosAng[2])
          + (2 * MIx[0] * MIx[2]) * (R_Len[0] * R_Len[2] * R_cosAng[1])
          + (2 * MIx[1] * MIx[2]) * (R_Len[1] * R_Len[2] * R_cosAng[0]);
      }
      //! d-spacing measure s = 1 / d = 2*sin(theta)/lambda.
      inline double s(const Miller::Index& MIx) const { return Q_as_s(Q(MIx));}
      //! d-spacing measure d = 1 / s = lamda/(2*sin(theta))
      inline double d(const Miller::Index& MIx) const { return Q_as_d(Q(MIx));}
      //@}

      //! @name Stream I/O.
      //@{
      //! Print the unit cell parameters to an output stream.
      std::ostream& print(std::ostream& os) const;
      //@}

    private:
      void SetVolume();
      void SetReciprocal();
      void SetMetricalMatrices();
      void SetOrthAndFracMatrix();
      void Initialize();

      Vec3   Len;
      Vec3   Ang;
      Vec3   sinAng;
      Vec3   cosAng;
      double Vol;
      Mx33   G;
      Vec3   R_Len;
      Vec3   R_Ang;
      Vec3   R_sinAng;
      Vec3   R_cosAng;
      Mx33   R_G;
      Mx33   Frac;
      Mx33   Orth;
  };

  //! iostream output operator for class UnitCell.
  std::ostream& operator<<(std::ostream& os, const UnitCell& uc);

} // namespace uctbx

#endif // CCTBX_UCTBX_H
