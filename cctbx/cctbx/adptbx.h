// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 11: Created (R.W. Grosse-Kunstleve)
 */

/*! \file
    Toolbox for the handling of anisotropic displacement parameters (ADPs).
 */

#ifndef CCTBX_ADPTBX_H
#define CCTBX_ADPTBX_H

#include <complex>
#include <utility>
#include <cctbx/uctbx.h>

namespace cctbx {
  //! ADP (anisotropic displacement parameters) Toolbox namespace.
  namespace adptbx {

  using namespace cctbx;

  static const error
    not_positive_definite(
      "anisotropic displacement tensor is not positive definite.");

  //! Constant: 2pi^2
  const double   TwoPiSquared = 2. * constants::pi * constants::pi;
  //! Constant: 8pi^2
  const double EightPiSquared = 8. * constants::pi * constants::pi;

  //! Convert isotropic displacement parameter U -> B.
  inline double
  U_as_B(double Uiso) {
    return Uiso * EightPiSquared;
  }
  //! Convert isotropic displacement parameter B -> U.
  inline double
  B_as_U(double Biso) {
    return Biso / EightPiSquared;
  }
  //! Convert anisotropic displacement parameters U -> B.
  template <class FloatType>
  carray<FloatType, 6>
  U_as_B(const carray<FloatType, 6>& Uaniso) {
    return EightPiSquared * Uaniso;
  }
  //! Convert anisotropic displacement parameters B -> U.
  template <class FloatType>
  carray<FloatType, 6>
  B_as_U(const carray<FloatType, 6>& Baniso) {
    return (1. / EightPiSquared) * Baniso;
  }

  //! Convert anisotropic displacement parameters Uuvrs -> Ustar.
  /*! The transformation matrix used is:<pre>
              (a*  0  0)
          C = ( 0 b*  0)
              ( 0  0 c*)</pre>
      The formula for the transformation is Ustar = C Uuvrs Ct,
      where Ct is the transposed of C. In this particular case
      the expression simplifies to:<pre>
          Ustar11 = a*^2  Uuvrs11
          Ustar22 = b*^2  Uuvrs22
          Ustar33 = c*^2  Uuvrs33
          Ustar12 = a* b* Uuvrs12
          Ustar13 = a* c* Uuvrs13
          Ustar23 = b* c* Uuvrs23</pre>
   */
  template <class FloatType>
  carray<FloatType, 6>
  Uuvrs_as_Ustar(const uctbx::UnitCell& uc,
                 const carray<FloatType, 6>& Uuvrs) {
    const double3& R_Len = uc.getLen(true);
    carray<FloatType, 6> Ustar;
    Ustar[0] = Uuvrs[0] * (R_Len[0] * R_Len[0]);
    Ustar[1] = Uuvrs[1] * (R_Len[1] * R_Len[1]);
    Ustar[2] = Uuvrs[2] * (R_Len[2] * R_Len[2]);
    Ustar[3] = Uuvrs[3] * (R_Len[0] * R_Len[1]);
    Ustar[4] = Uuvrs[4] * (R_Len[0] * R_Len[2]);
    Ustar[5] = Uuvrs[5] * (R_Len[1] * R_Len[2]);
    return Ustar;
  }
  //! Convert anisotropic displacement parameters Ustar -> Uuvrs.
  /*! Inverse of Uuvrs_as_Ustar().
   */
  template <class FloatType>
  carray<FloatType, 6>
  Ustar_as_Uuvrs(const uctbx::UnitCell& uc,
                 const carray<FloatType, 6>& Ustar) {
    const double3& R_Len = uc.getLen(true);
    carray<FloatType, 6> Uuvrs;
    Uuvrs[0] = Ustar[0] / (R_Len[0] * R_Len[0]);
    Uuvrs[1] = Ustar[1] / (R_Len[1] * R_Len[1]);
    Uuvrs[2] = Ustar[2] / (R_Len[2] * R_Len[2]);
    Uuvrs[3] = Ustar[3] / (R_Len[0] * R_Len[1]);
    Uuvrs[4] = Ustar[4] / (R_Len[0] * R_Len[2]);
    Uuvrs[5] = Ustar[5] / (R_Len[1] * R_Len[2]);
    return Uuvrs;
  }

  //! Convert anisotropic displacement parameters Ucart -> Ustar.
  /*! The transformation matrix C used is the orthogonalization
      matrix for the given UnitCell.
      The formula for the transformation is Ustar = C Ucart Ct,
      where Ct is the transposed of C.
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Ucart_as_Ustar(const uctbx::UnitCell& uc,
                 const carray<FloatType, 6>& Ucart) {
    return MatrixLite::CondensedTensorTransformation(
      uc.getFractionalizationMatrix(), Ucart);
  }
  //! Convert anisotropic displacement parameters Ustar -> Ucart.
  /*! Inverse of Ucart_as_Ustar(). I.e., the transformation matrix
      C used is the fractionalization matrix for the given
      UnitCell.
      The formula for the transformation is Ucart = C Ustar Ct,
      where Ct is the transposed of C.
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Ustar_as_Ucart(const uctbx::UnitCell& uc,
                 const carray<FloatType, 6>& Ustar) {
    return MatrixLite::CondensedTensorTransformation(
      uc.getOrthogonalizationMatrix(), Ustar);
  }

  //! Convert anisotropic displacement parameters Ucart -> Uuvrs.
  /*! This is implemented without a significant loss of efficiency
      as Ustar_as_Uuvrs(uc, Ucart_as_Ustar(uc, Ucart)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Ucart_as_Uuvrs(const uctbx::UnitCell& uc,
                 const carray<FloatType, 6>& Ucart) {
    return Ustar_as_Uuvrs(uc, Ucart_as_Ustar(uc, Ucart));
  }
  //! Convert anisotropic displacement parameters Uuvrs -> Ucart.
  /*! This is implemented without a significant loss of efficiency
      as Ustar_as_Ucart(uc, Uuvrs_as_Ustar(uc, Uuvrs)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Uuvrs_as_Ucart(const uctbx::UnitCell& uc,
                 const carray<FloatType, 6>& Uuvrs) {
    return Ustar_as_Ucart(uc, Uuvrs_as_Ustar(uc, Uuvrs));
  }

  //! Convert anisotropic displacement parameters Ustar -> beta.
  /*! The elements of Ustar are multiplied by 2pi^2.
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Ustar_as_beta(const carray<FloatType, 6>& Ustar) {
    return TwoPiSquared * Ustar;
  }
  //! Convert anisotropic displacement parameters beta -> Ustar.
  /*! The elements of beta are divided by 2pi^2.
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  beta_as_Ustar(const carray<FloatType, 6>& beta) {
    return beta / TwoPiSquared;
  }

  //! Convert anisotropic displacement parameters Ucart -> beta.
  /*! This is implemented as Ustar_as_beta(Ucart_as_Ustar(uc, Ucart)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Ucart_as_beta(const uctbx::UnitCell& uc,
                const carray<FloatType, 6>& Ucart) {
    return Ustar_as_beta(Ucart_as_Ustar(uc, Ucart));
  }
  //! Convert anisotropic displacement parameters beta -> Ucart.
  /*! This is implemented as Ustar_as_Ucart(uc, beta_as_Ustar(beta)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  beta_as_Ucart(const uctbx::UnitCell& uc,
                const carray<FloatType, 6>& beta) {
    return Ustar_as_Ucart(uc, beta_as_Ustar(beta));
  }

  //! Convert anisotropic displacement parameters Uuvrs -> beta.
  /*! This is implemented as Ustar_as_beta(Uuvrs_as_Ustar(uc, Uuvrs)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Uuvrs_as_beta(const uctbx::UnitCell& uc,
                const carray<FloatType, 6>& Uuvrs) {
    return Ustar_as_beta(Uuvrs_as_Ustar(uc, Uuvrs));
  }
  //! Convert anisotropic displacement parameters beta -> Uuvrs.
  /*! This is implemented as Ustar_as_Uuvrs(uc, beta_as_Ustar(beta)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  beta_as_Uuvrs(const uctbx::UnitCell& uc,
                const carray<FloatType, 6>& beta) {
    return Ustar_as_Uuvrs(uc, beta_as_Ustar(beta));
  }

  //! Convert anisotropic displacement parameters Ucart -> Uiso.
  /*! Uiso is defined as the mean of the diagonal elements of Ucart:<pre>
          Uiso = 1/3 (Ucart11 + Ucart22 + Ucart33)</pre>
   */
  template <class FloatType>
  inline FloatType
  Ucart_as_Uiso(const carray<FloatType, 6>& Ucart)
  {
    return (Ucart[0] + Ucart[1] + Ucart[2]) / 3.;
  }
  //! Convert Uiso -> anisotropic displacement parameters Ucart.
  /*! The diagonal elements of Ucart are set to the value of Uiso.
      The off-diagonal components Ucart are set to zero.
   */
  template <class FloatType>
  carray<FloatType, 6>
  Uiso_as_Ucart(const FloatType& Uiso)
  {
    carray<FloatType, 6> result;
    result.assign(0.);
    for(std::size_t i=0;i<3;i++) result[i] = Uiso;
    return result;
  }

  //! Convert Ustar -> Uiso.
  /*! This is implemented as Ucart_as_Uiso(Ustar_as_Ucart(uc, Ustar)).
   */
  template <class FloatType>
  inline FloatType
  Ustar_as_Uiso(const uctbx::UnitCell& uc,
                const carray<FloatType, 6>& Ustar)
  {
    return Ucart_as_Uiso(Ustar_as_Ucart(uc, Ustar));
  }
  //! Convert Uiso -> Ustar.
  /*! This is implemented as Ucart_as_Ustar(uc, Uiso_as_Ucart(Uiso)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Uiso_as_Ustar(const uctbx::UnitCell& uc,
                const FloatType& Uiso)
  {
    return Ucart_as_Ustar(uc, Uiso_as_Ucart(Uiso));
  }

  //! Convert Uuvrs -> Uiso.
  /*! This is implemented as Ucart_as_Uiso(Uuvrs_as_Ucart(uc, Uuvrs)).
   */
  template <class FloatType>
  inline FloatType
  Uuvrs_as_Uiso(const uctbx::UnitCell& uc,
                const carray<FloatType, 6>& Uuvrs)
  {
    return Ucart_as_Uiso(Uuvrs_as_Ucart(uc, Uuvrs));
  }
  //! Convert Uiso -> Uuvrs.
  /*! This is implemented as Ucart_as_Uuvrs(uc, Uiso_as_Ucart(Uiso)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Uiso_as_Uuvrs(const uctbx::UnitCell& uc,
                const FloatType& Uiso)
  {
    return Ucart_as_Uuvrs(uc, Uiso_as_Ucart(Uiso));
  }

  //! Convert beta -> Uiso.
  /*! This is implemented as Ucart_as_Uiso(beta_as_Ucart(uc, beta)).
   */
  template <class FloatType>
  inline FloatType
  beta_as_Uiso(const uctbx::UnitCell& uc,
               const carray<FloatType, 6>& beta)
  {
    return Ucart_as_Uiso(beta_as_Ucart(uc, beta));
  }
  //! Convert Uiso -> beta.
  /*! This is implemented as Ucart_as_beta(uc, Uiso_as_Ucart(Uiso)).
   */
  template <class FloatType>
  inline carray<FloatType, 6>
  Uiso_as_beta(const uctbx::UnitCell& uc,
                const FloatType& Uiso)
  {
    return Ucart_as_beta(uc, Uiso_as_Ucart(Uiso));
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
                         const carray<FloatType, 6>& Ustar)
  {
    return std::exp(-TwoPiSquared * (
        (MIx[0] * MIx[0]) * Ustar[0]
      + (MIx[1] * MIx[1]) * Ustar[1]
      + (MIx[2] * MIx[2]) * Ustar[2]
      + (2 * MIx[0] * MIx[1]) * Ustar[3]
      + (2 * MIx[0] * MIx[2]) * Ustar[4]
      + (2 * MIx[1] * MIx[2]) * Ustar[5]));
  }

  //! Anisotropic Debye-Waller factor given Miller index and beta.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactor_beta(const Miller::Index& MIx,
                         const carray<FloatType, 6>& beta)
  {
    return DebyeWallerFactorUstar(MIx, beta_as_Ustar(beta));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Uuvrs.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUuvrs(const uctbx::UnitCell& uc,
                         const Miller::Index& MIx,
                         const carray<FloatType, 6>& Uuvrs)
  {
    return DebyeWallerFactorUstar(MIx, Uuvrs_as_Ustar(uc, Uuvrs));
  }
  //! Anisotropic Debye-Waller factor given Miller index and Ucart.
  template <class FloatType>
  inline FloatType
  DebyeWallerFactorUcart(const uctbx::UnitCell& uc,
                         const Miller::Index& MIx,
                         const carray<FloatType, 6>& Ucart)
  {
    return DebyeWallerFactorUstar(MIx, Ucart_as_Ustar(uc, Ucart));
  }

  //! Determine the eigenvalues of the anisotropic displacement tensor.
  /*! Since the anisotropic displacement tensor is a symmetric matrix,
      all eigenvalues are real. The eigenvalues lambda are determined
      as the three real roots of the cubic equation
      <p>
          |(adp - lambda * I)| = 0
      <p>
      where I is the identity matrix. The solutions are
      obtained analytically using Cardan's formula.
      Detailed comments are embedded in the source code.
      <p>
      See also: Eigenvectors().
   */
  template <class FloatType>
  carray<FloatType, 3>
  Eigenvalues(const carray<FloatType, 6>& adp)
  {
    /* The eigenvalues lambda are found by:
         1. Determining the elements of the matrix (adp - lambda * I),
            where I is the identity matrix.
         2. Computing the determinant of that matrix as a function
            of lambda. This results in a cubic equation.
         3. Finding the real roots of the cubic equation with
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
    // to circumvent numerical instabilities due to rounding errors the
    // roots are determined as complex numbers.
    std::complex<FloatType> D(p * p * p / 27. + q * q / 4.);
    std::complex<FloatType> sqrtD = std::sqrt(D);
    FloatType mq2 = -q / 2.;
    std::complex<FloatType> u = std::pow(mq2 + sqrtD, 1/3.);
    std::complex<FloatType> v = std::pow(mq2 - sqrtD, 1/3.);
    std::complex<FloatType> epsilon1(-0.5, std::sqrt(3.) * 0.5);
    std::complex<FloatType> epsilon2 = std::conj(epsilon1);
    // since the anisotropic displacement tensor is a symmetric matrix,
    // all the imaginary components of the roots must be zero.
    carray<FloatType, 3> result;
    result[0] = (u + v).real();
    result[1] = (epsilon1 * u + epsilon2 * v).real();
    result[2] = (epsilon2 * u + epsilon1 * v).real();
    // convert the solutions y of the reduced form to the
    // solutions x of the normal form.
    for(std::size_t i = 0;i<3;i++) result[i] -= r / 3.;
    return result;
  }

  /*! \brief Test if the anisotropic displacement tensor is
      positive definite, given eigenvalues.
   */
  /*! Test if all eigenvalues are > 0.
      <p>
      See also: CheckPositiveDefinite(), Eigenvalues().
   */
  template <class FloatType>
  bool
  isPositiveDefinite(const carray<FloatType, 3>& adp_eigenvalues) {
    return adp_eigenvalues[carray_min_index(adp_eigenvalues)] > 0.;
  }

  /*! \brief Test if the anisotropic displacement tensor is
      positive definite.
   */
  /*! Test if all eigenvalues are > 0.
      <p>
      See also: CheckPositiveDefinite(), Eigenvalues().
   */
  template <class FloatType>
  bool
  isPositiveDefinite(const carray<FloatType, 6>& adp) {
    return isPositiveDefinite(Eigenvalues(adp));
  }

  /*! \brief Assert that the anisotropic displacement tensor is
      positive definite, given eigenvalues.
   */
  /*! An exception is thrown if the assertion fails.
      <p>
      See also: isPositiveDefinite(), Eigenvalues().
   */
  template <class FloatType>
  void
  CheckPositiveDefinite(const carray<FloatType, 3>& adp_eigenvalues) {
    if (!(isPositiveDefinite(adp_eigenvalues))) {
     throw not_positive_definite;
    }
  }

  /*! \brief Assert that the anisotropic displacement tensor is
      positive definite.
   */
  /*! An exception is thrown if the assertion fails.
      <p>
      See also: isPositiveDefinite(), Eigenvalues().
   */
  template <class FloatType>
  void
  CheckPositiveDefinite(const carray<FloatType, 6>& adp) {
    CheckPositiveDefinite(Eigenvalues(adp));
  }

  namespace detail {

    template <class FloatType>
    std::pair<carray<FloatType, 3>, FloatType>
    recursively_multiply(const carray<FloatType, 9>& M,
                         carray<FloatType, 3> V,
                         FloatType tolerance = 1.e-6)
    {
      unsigned int RunAwayCounter = 0;
      for (;;) {
        carray<FloatType, 3> MV;
        MatrixLite::multiply<FloatType>(M.elems, V.elems, 3, 3, 1, MV.elems);
        FloatType abs_lambda = std::sqrt(MV * MV);
        if (abs_lambda == 0.) throw not_positive_definite;
        MV = MV / abs_lambda;
        carray<FloatType, 3> absMV = carray_abs(MV);
        std::size_t iMax = carray_max_index(absMV);
        FloatType scaled_tolerance = absMV[iMax] * tolerance;
        bool converged = MatrixLite::approx_equal_scaled(
          MV, V, scaled_tolerance);
        if (!converged && MatrixLite::approx_equal_scaled(
          MV, -V, scaled_tolerance)) {
          throw not_positive_definite; // lambda < 0
        }
        V = MV;
        if (converged) return std::make_pair(V, abs_lambda);
        RunAwayCounter++;
        if (RunAwayCounter > 10000000) throw cctbx_internal_error();
      }
    }

  } // namespace detail

  //! Group of associated eigenvectors and values.
  template <class FloatType>
  class Eigensystem
  {
    public:
      //! Default constructor. Some data members are not initialized!
      Eigensystem() {}
      /*! \brief Determine the eigenvectors and eigenvalues of the
          anisotropic displacement tensor.
       */
      /*! Since the anisotropic displacement tensor is a symmetric matrix,
          all eigenvalues are real and the eigenvectors can be chosen
          orthonormal.
          <p>
          The eigenvectors are determined with the method of
          successive approximations as outlined in J.F. Nye,
          Physical Properties of Crystals, Oxford Science
          Publications, 1992, pp.165-168.
          The given tolerance is used to determine convergence.
          <p>
          An exception is thrown if any of the eigenvalues is
          less than or equal to zero. This indicates that the
          anisotropic displacement tensor is not positive definite.
          <p>
          See also: Eigenvalues().
       */
      Eigensystem(const carray<FloatType, 6>& adp,
                  FloatType tolerance = 1.e-6)
      {
        carray<FloatType, 9> M[2];
        M[0] = MatrixLite::CondensedSymMx33_as_FullSymMx33(adp,
               type_holder<FloatType>());
        FloatType d = MatrixLite::Determinant(M[0]);
        if (d == 0.) {
          throw not_positive_definite;
        }
        M[1] = MatrixLite::CoFactorMxTp(M[0]) / d;
        std::size_t iLarge[2];
        for(std::size_t iM=0;iM<2;iM++) {
          carray<FloatType, 3>
          absDiag = carray_abs(MatrixLite::DiagonalElements(M[iM]));
          iLarge[iM] = carray_max_index(absDiag);
          if (iM != 0 && iLarge[1] == iLarge[0]) {
            absDiag[iLarge[1]] = -1.;
            iLarge[1] = carray_max_index(absDiag);
            cctbx_assert(iLarge[1] != iLarge[0]);
          }
          carray<FloatType, 3> V;
          V.assign(0.);
          V[iLarge[iM]] = 1.;
          std::pair<carray<FloatType, 3>, FloatType>
          V_lambda = detail::recursively_multiply(M[iM], V);
          m_vectors[iM] = V_lambda.first;
          m_values[iM] = V_lambda.second;
        }
        m_vectors[2] = MatrixLite::cross_product(m_vectors[0], m_vectors[1]);
        cctbx_assert(m_vectors[2] * m_vectors[2] != 0.);
        m_values[1] = 1. / m_values[1];
        m_values[2] = (adp[0] + adp[1] + adp[2]) - (m_values[0] + m_values[1]);
      }
      //! Access the i'th eigenvector.
      /*! An exception is thrown if i >= 3.
       */
      const carray<FloatType, 3>&
      vectors(std::size_t i) const {
        if (i >= m_vectors.size()) throw error_index();
        return m_vectors[i];
      }
      //! Access the i'th eigenvalue.
      /*! An exception is thrown if i >= 3.
       */
      FloatType
      values(std::size_t i) const {
        if (i >= m_values.size()) throw error_index();
        return m_values[i];
      }
    private:
      carray<carray<FloatType, 3>, 3> m_vectors;
      carray<FloatType, 3> m_values;
  };

}} // namespace cctbx::adptbx

#endif // CCTBX_ADPTBX_H
