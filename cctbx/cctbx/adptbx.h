// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 11: Created (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for the handling of anisotropic displacement parameters.
 */

#ifndef CCTBX_ADPTBX_H
#define CCTBX_ADPTBX_H

#include <cctbx/uctbx.h>

//! ADP Toolbox namespace.
namespace adptbx {

  using namespace cctbx;

  static const error
    not_positive_definite("adp tensor is not positive definite.");

  const double   TwoPiSquared = 2. * constants::pi * constants::pi;
  const double EightPiSquared = 8. * constants::pi * constants::pi;

  template <class FloatType>
  struct return_type {};

  template <class FloatType6, class FloatType33>
  boost::array<FloatType33, 3*3>
  Xaniso_as_SymMx33(const boost::array<FloatType6, 6>& Xaniso,
                    return_type<FloatType33> t)
  {
    boost::array<FloatType33, 3*3> M;
    M[0] = Xaniso[0];
    M[1] = Xaniso[3];
    M[2] = Xaniso[4];
    M[3] = Xaniso[3];
    M[4] = Xaniso[1];
    M[5] = Xaniso[5];
    M[6] = Xaniso[4];
    M[7] = Xaniso[5];
    M[8] = Xaniso[2];
    return M;
  }

  template <class FloatType33, class FloatType6>
  inline boost::array<FloatType6, 6>
  SymMx33_as_Xaniso(const boost::array<FloatType33, 3*3>& M,
                    return_type<FloatType6> t)
  {
    boost::array<FloatType6, 6> Xaniso;
    Xaniso[0] = M[0];
    Xaniso[1] = M[4];
    Xaniso[2] = M[8];
    Xaniso[3] = M[1];
    Xaniso[4] = M[2];
    Xaniso[5] = M[5];
    return Xaniso;
  }

  template <class FloatType>
  boost::array<FloatType, 9>
  A_X_At(const boost::array<FloatType, 9>& A,
         const boost::array<FloatType, 9>& X)
  {
    boost::array<FloatType, 9> AX;
    MatrixLite::multiply<FloatType>(A.elems, X.elems, 3, 3, 3, AX.elems);
    boost::array<FloatType, 9> At;
    MatrixLite::transpose<FloatType>(A.elems, 3, 3, At.elems);
    boost::array<FloatType, 9> AXAt;
    MatrixLite::multiply<FloatType>(AX.elems, At.elems, 3, 3, 3, AXAt.elems);
    return AXAt;
  }

  template <class FloatType>
  inline boost::array<FloatType, 6>
  A_Xaniso_At(const uctbx::Mx33& A, const boost::array<FloatType, 6>& Xaniso)
  {
    boost::array<FloatType, 9>
    X = Xaniso_as_SymMx33(Xaniso, return_type<FloatType>());
    return SymMx33_as_Xaniso(A_X_At(A, X), return_type<FloatType>());
  }

  //! Convert isotropic adp U -> B.
  inline double
  U_as_B(double Uiso) {
    return Uiso * EightPiSquared;
  }
  //! Convert isotropic adp B -> U.
  inline double
  B_as_U(double Biso) {
    return Biso / EightPiSquared;
  }
  //! Convert anisotropic adp U -> B.
  template <class FloatType>
  boost::array<FloatType, 6>
  U_as_B(const boost::array<FloatType, 6>& Uaniso) {
    return EightPiSquared * Uaniso;
  }
  //! Convert anisotropic adp B -> U.
  template <class FloatType>
  boost::array<FloatType, 6>
  B_as_U(const boost::array<FloatType, 6>& Baniso) {
    return (1. / EightPiSquared) * Baniso;
  }

  //! Convert anisotropic adp Uuvrs -> Ustar.
  template <class FloatType>
  boost::array<FloatType, 6>
  Uuvrs_as_Ustar(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Uuvrs) {
    const uctbx::Vec3& R_Len = uc.getLen(true);
    boost::array<FloatType, 6> Ustar;
    Ustar[0] = Uuvrs[0] * (R_Len[0] * R_Len[0]);
    Ustar[1] = Uuvrs[1] * (R_Len[1] * R_Len[1]);
    Ustar[2] = Uuvrs[2] * (R_Len[2] * R_Len[2]);
    Ustar[3] = Uuvrs[3] * (R_Len[0] * R_Len[1]);
    Ustar[4] = Uuvrs[4] * (R_Len[0] * R_Len[2]);
    Ustar[5] = Uuvrs[5] * (R_Len[1] * R_Len[2]);
    return Ustar;
  }
  //! Convert anisotropic adp Ustar -> Uuvrs.
  template <class FloatType>
  boost::array<FloatType, 6>
  Ustar_as_Uuvrs(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Ustar) {
    const uctbx::Vec3& R_Len = uc.getLen(true);
    boost::array<FloatType, 6> Uuvrs;
    Uuvrs[0] = Ustar[0] / (R_Len[0] * R_Len[0]);
    Uuvrs[1] = Ustar[1] / (R_Len[1] * R_Len[1]);
    Uuvrs[2] = Ustar[2] / (R_Len[2] * R_Len[2]);
    Uuvrs[3] = Ustar[3] / (R_Len[0] * R_Len[1]);
    Uuvrs[4] = Ustar[4] / (R_Len[0] * R_Len[2]);
    Uuvrs[5] = Ustar[5] / (R_Len[1] * R_Len[2]);
    return Uuvrs;
  }

  //! Convert anisotropic adp Ucart -> Ustar.
  template <class FloatType>
  boost::array<FloatType, 6>
  Ucart_as_Ustar(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Ucart) {
    return A_Xaniso_At(uc.getFractionalizationMatrix(), Ucart);
  }
  //! Convert anisotropic adp Ustar -> Ucart.
  template <class FloatType>
  boost::array<FloatType, 6>
  Ustar_as_Ucart(const uctbx::UnitCell& uc,
                 const boost::array<FloatType, 6>& Ustar) {
    return A_Xaniso_At(uc.getOrthogonalizationMatrix(), Ustar);
  }

  //! Convert Uuvrs -> Uiso.
  // From Xtal 3.7.1 source code.
  template <class FloatType>
  inline FloatType
  Uuvrs_as_Uiso(const uctbx::UnitCell& uc,
                const boost::array<FloatType, 6>& Uuvrs)
  {
    const uctbx::Vec3&   Len = uc.getLen(false);
    const uctbx::Vec3& R_Len = uc.getLen(true);
    const uctbx::Vec3& cosAng = uc.get_cosAng(false);
    FloatType Uiso = 0.;
    FloatType LRL[3];
    for(std::size_t i=0;i<3;i++) {
      LRL[i] = Len[i] * R_Len[i];
      Uiso += Uuvrs[i] * LRL[i] * LRL[i];
    }
    Uiso += Uuvrs[3] * 2. * LRL[0] * LRL[1] * cosAng[2];
    Uiso += Uuvrs[4] * 2. * LRL[0] * LRL[2] * cosAng[1];
    Uiso += Uuvrs[5] * 2. * LRL[1] * LRL[2] * cosAng[0];
    return Uiso / 3.;
  }
  //! Convert Uiso -> Uuvrs.
  // From Xtal 3.7.1 source code.
  template <class FloatType>
  inline boost::array<FloatType, 6>
  Uiso_as_Uuvrs(const uctbx::UnitCell& uc,
                const FloatType& Uiso)
  {
    const uctbx::Vec3& R_cosAng = uc.get_cosAng(true);
    boost::array<FloatType, 6> Uuvrs;
    Uuvrs.assign(Uiso);
    Uuvrs[3] *= R_cosAng[2];
    Uuvrs[4] *= R_cosAng[1];
    Uuvrs[5] *= R_cosAng[0];
    return Uuvrs;
  }

  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and Biso.
  inline double
  DebyeWallerFactorBiso(double stol2,
                        double Biso)
  {
    return std::exp(-Biso * stol2);
  }
  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and Uiso.
  inline double
  DebyeWallerFactorUiso(double stol2,
                        double Uiso)
  {
    return DebyeWallerFactorBiso(stol2, U_as_B(Uiso));
  }
  //! Isotropic Debye-Waller factor given Miller index and Biso.
  inline double
  DebyeWallerFactorBiso(const uctbx::UnitCell& uc,
                        const Miller::Index& MIx,
                        double Biso)
  {
    return DebyeWallerFactorBiso(uc.Q(MIx) / 4., Biso);
  }
  //! Isotropic Debye-Waller factor given Miller index and Uiso.
  inline double
  DebyeWallerFactorUiso(const uctbx::UnitCell& uc,
                        const Miller::Index& MIx,
                        double Uiso)
  {
    return DebyeWallerFactorBiso(uc, MIx, U_as_B(Uiso));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Ustar.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUstar(const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Ustar)
  {
    return std::exp(-TwoPiSquared * (
        (MIx[0] * MIx[0]) * Ustar[0]
      + (MIx[1] * MIx[1]) * Ustar[1]
      + (MIx[2] * MIx[2]) * Ustar[2]
      + (2 * MIx[0] * MIx[1]) * Ustar[3]
      + (2 * MIx[0] * MIx[2]) * Ustar[4]
      + (2 * MIx[1] * MIx[2]) * Ustar[5]));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Uuvrs.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUuvrs(const uctbx::UnitCell& uc,
                         const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Uuvrs)
  {
    return DebyeWallerFactorUstar(MIx, Uuvrs_as_Ustar(uc, Uuvrs));
  }
  //! Anisotropic Debye-Waller factor given Miller index and Ucart.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUcart(const uctbx::UnitCell& uc,
                         const Miller::Index& MIx,
                         const boost::array<FloatType, 6>& Ucart)
  {
    return DebyeWallerFactorUstar(MIx, Ucart_as_Ustar(uc, Ucart));
  }

  //! Determine the eigenvalues of the adp tensor.
  /*! The eigenvalues lambda are determined as the three real roots
      of the cubic equation
      <p>
          |(adp - lambda * I)| = 0
      <p>
      where I is the identity matrix. The solution is
      obtained analytically using Cardan's formula.
      Detailed comments are embedded in the source code.
      <p>
      An exception is thrown if the cubic equation has imaginary
      roots. This indicates that the adp tensor is not
      positive definite.
   */
  template <class FloatType>
  boost::array<FloatType, 3>
  EigenValues(const boost::array<FloatType, 6>& adp)
  {
    /* The eigenvalues lambda are found by:
         1. Determining the elements of the matrix (adp - lambda * I),
            where I is the identity matrix.
         2. Computing the determinant of that matrix as a function
            of lambda. This results in a cubic equation.
         3. Finding the three real roots of the cubic equation with
            Cardan's formula, taken from Taschenbuch der Mathematik
            by Bronstein & Semendjajew (p. 183 in the reprint of
            the 20th edition from 1983).
    */
    /* Mathematica code used to generate the coefficients of the
       "normal form" of the cubic equation:
         U={{adp0,adp3,adp4},{adp3,adp1,adp5},{adp4,adp5,adp2}}
         L={{lambda,0,0},{0,lambda,0},{0,0,lambda}}
         Det[U-L]*(-1)
         CForm[Det[U-L]*(-1)]
     */
    // normal form: x^3 + r x^2 + s x + t == 0
    FloatType r = -adp[0] - adp[1] - adp[2];
    FloatType s =   adp[0] * adp[1] + adp[0] * adp[2] + adp[1] * adp[2]
                  - adp[3] * adp[3] - adp[4] * adp[4] - adp[5] * adp[5];
    FloatType t =   adp[0] * adp[5] * adp[5] - adp[0] * adp[1] * adp[2]
                  + adp[2] * adp[3] * adp[3] + adp[1] * adp[4] * adp[4]
                  - 2 * adp[3] * adp[4] * adp[5];
    /* "Reduced form" of the cubic equation according to
       Taschenbuch der Mathematik.
     */
    // reduced form: y^3 + p y + q == 0
    FloatType p = s - r * r / 3.;
    FloatType q = 2. * r * r * r / 27. - r * s / 3. + t;
    /* "D" is computed to test for the number of real roots.
       There are three real roots only if D < 0.
     */
    FloatType p3 = p * p * p;
    FloatType q2 = q * q;
    FloatType D = p3 / 27. + q2 / 4.;
    if (D >= 0.) throw not_positive_definite;
    /* At this point it is clear that there are three real eigenvaues.
       They can now be determined without involving complex numbers.
     */
    FloatType zeta = std::sqrt(-p3 / 27);
    FloatType phi = std::acos(-q / (2. * zeta));
    FloatType rzeta = 2. * std::pow(zeta, FloatType(1/3.));
    boost::array<FloatType, 3> result;
    for (int i = 0; i < 3; i++) {
      result[i] =   rzeta * std::cos((phi + (i * 2) * constants::pi) / 3.)
                  - r / 3.; // the term on this line converts the
    }                       // solutions y of the reduced form to
    return result;          // the solutions x of the normal form.
  }

  namespace detail {

    template <class FloatType>
    boost::array<FloatType, 3>
    recursively_multiply(const boost::array<FloatType, 9>& M,
                         boost::array<FloatType, 3> V,
                         FloatType tolerance = 1.e-6)
    {
      for (;;) {
        boost::array<FloatType, 3> MV;
        MatrixLite::multiply<FloatType>(M.elems, V.elems, 3, 3, 1, MV.elems);
        FloatType LenMV = std::sqrt(MV * MV);
        if (LenMV == 0.) return MV;
        MV = MV / LenMV;
        boost::array<FloatType, 3> absMV = boost::array_abs(MV);
        std::size_t iMax = boost::array_max_index(absMV);
        FloatType scaled_tolerance = absMV[iMax] * tolerance;
        bool converged = MatrixLite::approx_equal(MV, V, scaled_tolerance);
        V = MV;
        if (converged) break;
      }
      return V;
    }

  } // namespace detail

  //! Determine the eigenvalues of the adp tensor.
  /*! XXX
   */
  template <class FloatType>
  boost::array<boost::array<FloatType, 3>, 3>
  EigenVectors(const boost::array<FloatType, 6>& adp, double tolerance = 1.e-6)
  {
    boost::array<FloatType, 9> M[2];
    M[0] = Xaniso_as_SymMx33(adp, return_type<FloatType>());
    FloatType d = MatrixLite::Determinant(M[0]);
    if (d == 0.) {
      throw not_positive_definite;
    }
    M[1] = MatrixLite::CoFactorMxTp(M[0]) / d;
    boost::array<boost::array<FloatType, 3>, 3> result;
    for(std::size_t iM=0;iM<2;iM++) {
      boost::array<FloatType, 3>
      Diag = MatrixLite::DiagonalElements(M[iM]);
      std::size_t iLarge = boost::array_max_index(boost::array_abs(Diag));
      boost::array<FloatType, 3> V;
      V.assign(0.);
      V[iLarge] = 1.;
      result[iM] = detail::recursively_multiply(M[iM], V);
      if (result[iM] * result[iM] == 0.) {
        throw not_positive_definite;
      }
    }
    result[2] = MatrixLite::cross_product(result[0], result[1]);
    cctbx_assert(result[2] * result[2] != 0.);
    return result;
  }

} // namespace adptbx

#endif // CCTBX_ADPTBX_H
