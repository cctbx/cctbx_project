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

#if defined(__GNUC__) && __GNUC__ < 3
# include <iostream>
#else
# include <ostream>
#endif
#include <cctbx/fixes/cmath>
#include <cctbx/error.h>
#include <cctbx/basic/matrixlite.h>
#include <cctbx/coordinates.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/matrix.h>
#include <cctbx/constants.h>
#include <cctbx/array_family/shared.h>

namespace cctbx {
  //! Unit Cell Toolbox namespace.
  namespace uctbx {

  static const error
    corrupt_unit_cell_parameters("Corrupt unit cell parameters.");
  static const error
    corrupt_metrical_matrix("Corrupt metrical matrix.");

  //! @name Conversions between radians and degrees.
  //@{
  inline double deg_as_rad(double deg) { return deg * constants::pi_180;}
  inline double rad_as_deg(double rad) { return rad / constants::pi_180;}
  //@}

  //! inline function for fast matrix * vector computation.
  template <class FloatType>
  inline af::tiny<FloatType, 3>
  operator*(const af::double9& m, const af::tiny<FloatType, 3>& v)
  {
    af::tiny<FloatType, 3> mv;
    mv[0] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
    mv[1] = m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
    mv[2] = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
    return mv;
  }

  //! Helper function for transforming metrical matrices.
  inline af::double9 getRtGR(const af::double9& G, const af::double9& R)
  {
    af::double9 Rt, GR, RtGR;
    MatrixLite::transpose<double>(R.begin(), 3, 3, Rt.begin());
    MatrixLite::multiply<double>(G.begin(), R.begin(), 3,3,3, GR.begin());
    MatrixLite::multiply<double>(Rt.begin(), GR.begin(), 3,3,3, RtGR.begin());
    return RtGR;
  }

  //! @name inline functions for fast conversions of d-spacing measures.
  //@{
  inline double Q_as_stol2(double Q) { return Q * .25; }
  inline double Q_as_two_stol(double Q) { return std::sqrt(Q); }
  inline double Q_as_stol(double Q) { return std::sqrt(Q) * .5; }
  inline double Q_as_two_theta(double Q, double wavelength, bool deg = false)
  {
    double result = 2. * std::asin(Q_as_stol(Q) * wavelength);
    if (deg) return rad_as_deg(result);
    return result;
  }
  inline double Q_as_d(double Q) {
    if (Q == 0.) return -1.;
    /* else */   return 1. / std::sqrt(Q);
  }
  //@}

  //! Helper class for passing unit cell parameters.
  class uc_params : public af::double6 {
    public:
      //! @name Constructors.
      //@{
      //! Constructor using parameters (a, b, c, alpha, beta, gamma).
      explicit
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
      uc_params(const af::double3& Len, const af::double3& Ang) {
        int i;
        for(i=0;i<3;i++) elems[i]     = Len[i];
        for(i=0;i<3;i++) elems[i + 3] = Ang[i];
      }
      //@}
      //! @name Access to arrays of lengths and angles.
      //@{
      //!
      double* Len() { return &elems[0]; }
      double* Ang() { return &elems[3]; }
      const double* Len() const { return &elems[0]; }
      const double* Ang() const { return &elems[3]; }
      //@}
      //! @name Access to individual lengths and angles.
      //@{
      double Len(int i) const { return elems[i]; }
      double Ang(int i) const { return elems[3 + i]; }
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
      explicit
      UnitCell(const uc_params& ucp);
      //! Constructor using parameters derived from a metrical matrix.
      /*! The metrical matrix is defined as:
         <pre>
         ( a*a,            a*b*cos(gamma), a*c*cos(beta)  )
         ( a*b*cos(gamma), b*b,            b*c*cos(alpha) )
         ( a*c*cos(beta),  b*c*cos(alpha), c*c            )
         </pre>
       */
      explicit
      UnitCell(const af::double9& MetricalMatrix);
      //@}

      //! @name Query parameters and volume.
      //@{
      uc_params getParameters(bool reciprocal = false) const;
      const af::double3& getLen(bool reciprocal = false) const {
        if (reciprocal == false) return Len;
        /* else */               return R_Len;
      }
      const af::double3& getAng(bool reciprocal = false) const {
        if (reciprocal == false) return Ang;
        /* else */               return R_Ang;
      }
      const af::double3& get_sinAng(bool reciprocal = false) const {
        if (reciprocal == false) return sinAng;
        /* else */               return R_sinAng;
      }
      const af::double3& get_cosAng(bool reciprocal = false) const {
        if (reciprocal == false) return cosAng;
        /* else */               return R_cosAng;
      }
      const af::double9& getMetricalMatrix(bool reciprocal = false) const {
        if (reciprocal == false) return G;
        /* else */               return R_G;
      }
      double getVolume() const { return Vol; }
      //@}

      //! @name Length^2 of the longest lattice vector in the unit cell.
      //@{
      double getLongestVector2() const { return LongestVector2; }
      //@}

      //! @name Test equality.
      //@{
      //! Test the equality of two Unit Cell instances.
      /*! Test if
          2 * abs(  (self.paramter - other.paramter)
                  / (self.paramter + other.paramter))
          is less then the given tolerance for all six unit cell
          parameters.
       */
      bool isEqual(const UnitCell& other, double tolerance = 1.e-6) const;
      //@}

      //! @name Orthogonalization and fractionalization of coordinates.
      //@{
      //! This matrix converts cartesian to fractional coordinates.<br>
      //! x(fractional) = matrix * x(cartesian).
      const af::double9& getFractionalizationMatrix() const { return Frac; }
      //! This matrix converts fractional to cartesian coordinates.<br>
      //! x(cartesian) = matrix * x(fractional).
      const af::double9& getOrthogonalizationMatrix() const { return Orth; }
      //! Converts cartesian coordinates Xc to fractional coordinates.
      template <class FloatType>
      fractional<FloatType>
      fractionalize(const cartesian<FloatType>& Xc) const {
        return Frac * Xc;
      }
      //! Converts fractional coordinates Xf to cartesian coordinates.
      template <class FloatType>
      cartesian<FloatType>
      orthogonalize(const fractional<FloatType>& Xf) const {
        return Orth * Xf;
      }
      //@}

      //! @name Measurements, given fractional coordinates.
      //@{
      //! Length squared of vector.
      template <class FloatType>
      FloatType Length2(const fractional<FloatType>& Xf) const {
        return orthogonalize(Xf).Length2();
      }
      //! Length of vector.
      template <class FloatType>
      FloatType Length(const fractional<FloatType>& Xf) const {
        return std::sqrt(Length2(Xf));
      }
      //! Distance squared.
      template <class FloatType>
      FloatType Distance2(const fractional<FloatType>& Xf,
                          const fractional<FloatType>& Yf) const {
        return Length2(fractional<FloatType>(Xf - Yf));
      }
      //! Distance.
      template <class FloatType>
      FloatType Distance(const fractional<FloatType>& Xf,
                         const fractional<FloatType>& Yf) const {
        return Length(fractional<FloatType>(Xf - Yf));
      }
      //! Shortest length squared under applicaton of periodicity.
      template <class FloatType>
      FloatType modShortLength2(const af::tiny<FloatType, 3>& Xf) const {
        return Length2(fractional<FloatType>(Xf).modShort());
      }
      //! Shortest length under applicaton of periodicity.
      template <class FloatType>
      FloatType modShortLength(const af::tiny<FloatType, 3>& Xf) const {
        return std::sqrt(modShortLength2(Xf));
      }
      //! Shortest distance squared under applicaton of periodicity.
      template <class FloatType>
      FloatType modShortDistance2(const fractional<FloatType>& Xf,
                                  const fractional<FloatType>& Yf) const {
        return modShortLength2(fractional<FloatType>(Xf - Yf));
      }
      //! Shortest distance under applicaton of periodicity.
      template <class FloatType>
      FloatType modShortDistance(const fractional<FloatType>& Xf,
                                 const fractional<FloatType>& Yf) const {
        return std::sqrt(modShortDistance2(Xf, Yf));
      }
      //@}

      //! @name Transformation (change-of-basis) of unit cell parameters.
      //@{
      //! InvCBMxR is the inverse of the 3x3 change-of-basis matrix
      //! that transforms coordinates in the old basis system to
      //! coodinates in the new basis system.
      UnitCell ChangeBasis(const af::double9& InvCBMxR, double RBF = 1.) const;
      UnitCell ChangeBasis(const sgtbx::RotMx& InvCBMxR) const;
      //@}

      //! @name Methods using Miller indices.
      //@{
      //! Compute the maximum Miller indices for a given minimum d-spacing.
      miller::Index MaxMillerIndices(double dmin) const;
      //! d-spacing measure Q = 1/d^2 = (2*sin(theta)/lambda)^2.
      double
      Q(const miller::Index& MIx) const
      {
        return
            (MIx[0] * MIx[0]) * R_G[0]
          + (MIx[1] * MIx[1]) * R_G[4]
          + (MIx[2] * MIx[2]) * R_G[8]
          + (2 * MIx[0] * MIx[1]) * R_G[1]
          + (2 * MIx[0] * MIx[2]) * R_G[2]
          + (2 * MIx[1] * MIx[2]) * R_G[5];
      }
      //! d-spacing measure Q = 1/d^2 = (2*sin(theta)/lambda)^2.
      af::shared<double>
      Q(const af::shared<miller::Index>& MIx) const;
      //! Maximum Q for given list of Miller indices.
      double
      max_Q(const af::shared<miller::Index>& MIx) const;
      //! Minimum and maximum Q for given list of Miller indices.
      af::double2
      min_max_Q(const af::shared<miller::Index>& MIx) const;

      //! d-spacing measure (sin(theta)/lambda)^2 = Q/4.
      double
      stol2(const miller::Index& MIx) const {
        return Q_as_stol2(Q(MIx));
      }
      //! d-spacing measure (sin(theta)/lambda)^2 = Q/4.
      af::shared<double>
      stol2(const af::shared<miller::Index>& MIx) const;

      //! d-spacing measure 2*sin(theta)/lambda = 1/d = sqrt(Q).
      double
      two_stol(const miller::Index& MIx) const {
        return Q_as_two_stol(Q(MIx));
      }
      //! d-spacing measure 2*sin(theta)/lambda = 1/d = sqrt(Q).
      af::shared<double>
      two_stol(const af::shared<miller::Index>& MIx) const;

      //! d-spacing measure sin(theta)/lambda = 1/(2*d) = sqrt(Q)/2.
      double
      stol(const miller::Index& MIx) const {
        return Q_as_stol(Q(MIx));
      }
      //! d-spacing measure sin(theta)/lambda = 1/(2*d) = sqrt(Q)/2.
      af::shared<double>
      stol(const af::shared<miller::Index>& MIx) const;

      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      double
      d(const miller::Index& MIx) const { return Q_as_d(Q(MIx));}
      //! d-spacing measure d = 1/(2*sin(theta)/lambda).
      af::shared<double>
      d(const af::shared<miller::Index>& MIx) const;

      //! Diffraction angle 2-theta in degrees, given wavelength lamdda.
      double
      two_theta(
        const miller::Index& MIx, double wavelength, bool deg = false) const
      {
        return Q_as_two_theta(Q(MIx), wavelength);
      }
      //! Diffraction angle 2-theta in degrees, given wavelength lamdda.
      af::shared<double>
      two_theta(
        af::shared<miller::Index> MIx,
        double wavelength,
        bool deg = false) const;
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

      af::double3 Len;
      af::double3 Ang;
      af::double3 sinAng;
      af::double3 cosAng;
      double      Vol;
      af::double9 G;
      af::double3 R_Len;
      af::double3 R_Ang;
      af::double3 R_sinAng;
      af::double3 R_cosAng;
      af::double9 R_G;
      af::double9 Frac;
      af::double9 Orth;
      double      LongestVector2;
  };

  //! ostream output operator for class UnitCell.
  std::ostream& operator<<(std::ostream& os, const UnitCell& uc);

  //! XXX
  template <typename FloatType>
  class incremental_d_star_sq
  {
    public:
      incremental_d_star_sq() {}

      incremental_d_star_sq(const UnitCell& ucell)
      {
        initialize(ucell.getMetricalMatrix(true));
      }

      void update0(int h0)
      {
        h0_ = h0;
        im0_ = (h0_ * h0_) * r_g00_;
      }

      void update1(int h1)
      {
        h1_ = h1;
        im1_ = im0_ + (h1_ * h1_) * r_g11_
                    + (2 * h0_ * h1_) * r_g01_;
      }

      FloatType get(int h2)
      {
        return im1_ + (h2 * h2) * r_g22_
                    + (2 * h0_ * h2) * r_g02_
                    + (2 * h1_ * h2) * r_g12_;
      }

    protected:
      FloatType r_g00_, r_g11_, r_g22_, r_g01_, r_g02_, r_g12_;
      int h0_, h1_;
      FloatType im0_, im1_;

      void initialize(const af::tiny<FloatType, 9>& r_g)
      {
        r_g00_ = r_g[0];
        r_g11_ = r_g[4];
        r_g22_ = r_g[8];
        r_g01_ = r_g[1];
        r_g02_ = r_g[2];
        r_g12_ = r_g[5];
      }
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_H
