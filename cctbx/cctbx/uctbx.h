// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
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
#include <cctbx/coordinates.h>
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
  template <class T>
  inline boost::array<T, 3>
  operator*(const Mx33& m, const boost::array<T, 3>& v)
  {
    boost::array<T, 3> mv;
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
      //! Default (1, 1, 1, 90, 90, 90).
      UnitCell();
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

      //! @name Length^2 of the longest lattice vector in the unit cell.
      //@{
      inline double getLongestVector2() const { return LongestVector2; }
      //@}

      //! @name Test equality.
      //@{
      //! Test the equality of two Unit Cell instances.  Test the fractional
      //! difference of each of the six parameters and compare to epsilon
      bool isEqual(const UnitCell& uc, const double& epsilon) const;
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
      template <class T>
      inline fractional<T>
      fractionalize(const cartesian<T>& Xc) const {
        return Frac * Xc;
      }
      //! Converts fractional coordinates Xf to cartesian coordinates.
      template <class T>
      inline cartesian<T>
      orthogonalize(const fractional<T>& Xf) const {
        return Orth * Xf;
      }
      //@}

      //! @name Measurements, given fractional coordinates.
      //@{
      //! Length squared of vector.
      template <class T>
      inline T Length2(const fractional<T>& Xf) const {
        return orthogonalize(Xf).Length2();
      }
      //! Length of vector.
      template <class T>
      inline T Length(const fractional<T>& Xf) const {
        return std::sqrt(Length2(Xf));
      }
      //! Distance squared.
      template <class T>
      inline T Distance2(const fractional<T>& Xf,
                         const fractional<T>& Yf) const {
        return Length2(fractional<T>(Xf - Yf));
      }
      //! Distance.
      template <class T>
      inline T Distance(const fractional<T>& Xf,
                        const fractional<T>& Yf) const {
        return Length(fractional<T>(Xf - Yf));
      }
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
      //! Inverse operation of MaxMiller indices.
      double MaxResolution(Miller::Index) const;
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

      //! @name Temperature factors.
      /*! Literature:
          Jan Drenth,
          Principles of Protein X-Ray Crystallography,
          2nd edition, 1999, section 4.9 (pp. 89-92)
       */
      //@{
      //! Convert isotropic temperature factor coefficient Uiso -> Biso.
      inline double Uiso_as_Biso(double Uiso) const
      {
        using cctbx::constants::pi;
        return Uiso * (8. * pi * pi);
      }
      //! Convert isotropic temperature factor coefficient Biso -> Uiso.
      inline double Biso_as_Uiso(double Biso) const
      {
        using cctbx::constants::pi;
        return Biso / (8. * pi * pi);
      }
      //! Convert anisotropic temperature factor coefficients Uij -> Bij.
      template <class FloatType>
      boost::array<FloatType, 6>
      Uij_as_Bij(const boost::array<FloatType, 6>& Uij) const
      {
        using cctbx::constants::pi;
        const FloatType TpiS = 2. * pi * pi;
        boost::array<FloatType, 6> Bij;
        Bij[0] = Uij[0] * (TpiS * (R_Len[0] * R_Len[0]));
        Bij[1] = Uij[1] * (TpiS * (R_Len[1] * R_Len[1]));
        Bij[2] = Uij[2] * (TpiS * (R_Len[2] * R_Len[2]));
        Bij[3] = Uij[3] * (TpiS * (R_Len[0] * R_Len[1]));
        Bij[4] = Uij[4] * (TpiS * (R_Len[0] * R_Len[2]));
        Bij[5] = Uij[5] * (TpiS * (R_Len[1] * R_Len[2]));
        return Bij;
      }
      //! Convert anisotropic temperature factor coefficients Bij -> Uij.
      template <class FloatType>
      boost::array<FloatType, 6>
      Bij_as_Uij(const boost::array<FloatType, 6>& Bij) const
      {
        using cctbx::constants::pi;
        const FloatType TpiS = 2. * pi * pi;
        boost::array<FloatType, 6> Uij;
        Uij[0] = Bij[0] / (TpiS * (R_Len[0] * R_Len[0]));
        Uij[1] = Bij[1] / (TpiS * (R_Len[1] * R_Len[1]));
        Uij[2] = Bij[2] / (TpiS * (R_Len[2] * R_Len[2]));
        Uij[3] = Bij[3] / (TpiS * (R_Len[0] * R_Len[1]));
        Uij[4] = Bij[4] / (TpiS * (R_Len[0] * R_Len[2]));
        Uij[5] = Bij[5] / (TpiS * (R_Len[1] * R_Len[2]));
        return Uij;
      }
      //! Convert Uij -> Uiso.
      template <class FloatType>
      inline FloatType
      Uij_as_Uiso(const boost::array<FloatType, 6>& Uij) const
      {
        FloatType Uiso = 0.;
        FloatType LRL[3];
        for(std::size_t i=0;i<3;i++) {
          LRL[i] = Len[i] * R_Len[i];
          Uiso += Uij[i] * LRL[i] * LRL[i];
        }
        Uiso += Uij[3] * 2. * LRL[0] * LRL[1] * cosAng[2];
        Uiso += Uij[4] * 2. * LRL[0] * LRL[2] * cosAng[1];
        Uiso += Uij[5] * 2. * LRL[1] * LRL[2] * cosAng[0];
        return Uiso / 3.;
      }
      //! Convert Uiso -> Uij.
      template <class FloatType>
      inline boost::array<FloatType, 6>
      Uiso_as_Uij(const FloatType& Uiso) const
      {
        boost::array<FloatType, 6> Uij;
        Uij.assign(Uiso);
        Uij[3] *= R_cosAng[2];
        Uij[4] *= R_cosAng[1];
        Uij[5] *= R_cosAng[0];
        return Uij;
      }
      //! Isotropic temperature factor given (sin(theta)/lambda)^2 and Biso.
      inline double
      TemperatureFactorB(double stol2,
                         double Biso) const
      {
        using cctbx::constants::pi;
        return std::exp(-Biso * stol2);
      }
      //! Isotropic temperature factor given (sin(theta)/lambda)^2 and Uiso.
      inline double
      TemperatureFactorU(double stol2,
                         double Uiso) const
      {
        return TemperatureFactorB(stol2, Uiso_as_Biso(Uiso));
      }
      //! Isotropic temperature factor given Miller index and Biso.
      inline double
      TemperatureFactorB(const Miller::Index& MIx,
                         double Biso) const
      {
        return TemperatureFactorB(Q(MIx) / 4., Biso);
      }
      //! Isotropic temperature factor given Miller index and Uiso.
      inline double
      TemperatureFactorU(const Miller::Index& MIx,
                         double Uiso) const
      {
        return TemperatureFactorB(MIx, Uiso_as_Biso(Uiso));
      }
      //! Anisotropic temperature factor given coefficients Bij.
      template <class FloatType>
      inline FloatType
      TemperatureFactorB(const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Bij) const
      {
        using cctbx::constants::pi;
        return std::exp(-(
            (MIx[0] * MIx[0]) * Bij[0]
          + (MIx[1] * MIx[1]) * Bij[1]
          + (MIx[2] * MIx[2]) * Bij[2]
          + (2 * MIx[0] * MIx[1]) * Bij[3]
          + (2 * MIx[0] * MIx[2]) * Bij[4]
          + (2 * MIx[1] * MIx[2]) * Bij[5]));
      }
      //! Anisotropic temperature factor given coefficients Uij.
      template <class FloatType>
      inline FloatType
      TemperatureFactorU(const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Uij) const
      {
        using cctbx::constants::pi;
        return std::exp(FloatType(-2. * pi * pi) * (
            (MIx[0] * MIx[0]) * (R_Len[0] * R_Len[0]) * Uij[0]
          + (MIx[1] * MIx[1]) * (R_Len[1] * R_Len[1]) * Uij[1]
          + (MIx[2] * MIx[2]) * (R_Len[2] * R_Len[2]) * Uij[2]
          + (2 * MIx[0] * MIx[1]) * (R_Len[0] * R_Len[1]) * Uij[3]
          + (2 * MIx[0] * MIx[2]) * (R_Len[0] * R_Len[2]) * Uij[4]
          + (2 * MIx[1] * MIx[2]) * (R_Len[1] * R_Len[2]) * Uij[5]));
      }
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
      void SetLongestVector2();
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
      double LongestVector2;
  };

  //! iostream output operator for class UnitCell.
  std::ostream& operator<<(std::ostream& os, const UnitCell& uc);

} // namespace uctbx

#endif // CCTBX_UCTBX_H
