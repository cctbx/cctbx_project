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
#include <cctbx/array_family/tiny_algebra.h>
#include <cctbx/uctbx.h>

namespace cctbx {
  //! ADP (anisotropic displacement parameters) Toolbox namespace.
  namespace adptbx {

  using namespace cctbx;

  static const error
    not_positive_definite(
      "anisotropic displacement tensor is not positive definite.");

  //! Convert isotropic displacement parameter U -> B.
  inline double
  U_as_B(double Uiso) {
    return Uiso * constants::eight_pi_sq;
  }
  //! Convert isotropic displacement parameter B -> U.
  inline double
  B_as_U(double Biso) {
    return Biso / constants::eight_pi_sq;
  }
  //! Convert anisotropic displacement parameters U -> B.
  template <typename FloatType>
  af::tiny<FloatType, 6>
  U_as_B(const af::tiny_plain<FloatType, 6>& Uaniso) {
    return constants::eight_pi_sq * af::tiny<FloatType, 6>(Uaniso);
  }
  //! Convert anisotropic displacement parameters B -> U.
  template <typename FloatType>
  af::tiny<FloatType, 6>
  B_as_U(const af::tiny<FloatType, 6>& Baniso) {
    return (1. / constants::eight_pi_sq) * Baniso;
  }

  //! Convert anisotropic displacement parameters Ucif -> Ustar.
  /*! The transformation matrix used is:<pre>
              (a*  0  0)
          C = ( 0 b*  0)
              ( 0  0 c*)</pre>
      The formula for the transformation is Ustar = C Ucif Ct,
      where Ct is the transposed of C. In this particular case
      the expression simplifies to:<pre>
          Ustar11 = a*^2  Ucif11
          Ustar22 = b*^2  Ucif22
          Ustar33 = c*^2  Ucif33
          Ustar12 = a* b* Ucif12
          Ustar13 = a* c* Ucif13
          Ustar23 = b* c* Ucif23</pre>
   */
  template <typename FloatType>
  af::tiny<FloatType, 6>
  Ucif_as_Ustar(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& Ucif) {
    const af::double3& R_Len = uc.getLen(true);
    af::tiny<FloatType, 6> Ustar;
    Ustar[0] = Ucif[0] * (R_Len[0] * R_Len[0]);
    Ustar[1] = Ucif[1] * (R_Len[1] * R_Len[1]);
    Ustar[2] = Ucif[2] * (R_Len[2] * R_Len[2]);
    Ustar[3] = Ucif[3] * (R_Len[0] * R_Len[1]);
    Ustar[4] = Ucif[4] * (R_Len[0] * R_Len[2]);
    Ustar[5] = Ucif[5] * (R_Len[1] * R_Len[2]);
    return Ustar;
  }
  //! Convert anisotropic displacement parameters Ustar -> Ucif.
  /*! Inverse of Ucif_as_Ustar().
   */
  template <typename FloatType>
  af::tiny<FloatType, 6>
  Ustar_as_Ucif(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& Ustar) {
    const af::double3& R_Len = uc.getLen(true);
    af::tiny<FloatType, 6> Ucif;
    Ucif[0] = Ustar[0] / (R_Len[0] * R_Len[0]);
    Ucif[1] = Ustar[1] / (R_Len[1] * R_Len[1]);
    Ucif[2] = Ustar[2] / (R_Len[2] * R_Len[2]);
    Ucif[3] = Ustar[3] / (R_Len[0] * R_Len[1]);
    Ucif[4] = Ustar[4] / (R_Len[0] * R_Len[2]);
    Ucif[5] = Ustar[5] / (R_Len[1] * R_Len[2]);
    return Ucif;
  }

  //! Convert anisotropic displacement parameters Ucart -> Ustar.
  /*! The transformation matrix C used is the orthogonalization
      matrix for the given UnitCell.
      The formula for the transformation is Ustar = C Ucart Ct,
      where Ct is the transposed of C.
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ucart_as_Ustar(const uctbx::UnitCell& uc,
                 const af::tiny<FloatType, 6>& Ucart) {
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
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ustar_as_Ucart(const uctbx::UnitCell& uc,
                 const af::tiny<FloatType, 6>& Ustar) {
    return MatrixLite::CondensedTensorTransformation(
      uc.getOrthogonalizationMatrix(), Ustar);
  }

  //! Convert anisotropic displacement parameters Ucart -> Ucif.
  /*! This is implemented without a significant loss of efficiency
      as Ustar_as_Ucif(uc, Ucart_as_Ustar(uc, Ucart)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ucart_as_Ucif(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& Ucart) {
    return Ustar_as_Ucif(uc, Ucart_as_Ustar(uc, Ucart));
  }
  //! Convert anisotropic displacement parameters Ucif -> Ucart.
  /*! This is implemented without a significant loss of efficiency
      as Ustar_as_Ucart(uc, Ucif_as_Ustar(uc, Ucif)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ucif_as_Ucart(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& Ucif) {
    return Ustar_as_Ucart(uc, Ucif_as_Ustar(uc, Ucif));
  }

  //! Convert anisotropic displacement parameters Ustar -> beta.
  /*! The elements of Ustar are multiplied by 2pi^2.
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ustar_as_beta(const af::tiny<FloatType, 6>& Ustar) {
    return constants::two_pi_sq * Ustar;
  }
  //! Convert anisotropic displacement parameters beta -> Ustar.
  /*! The elements of beta are divided by 2pi^2.
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  beta_as_Ustar(const af::tiny<FloatType, 6>& beta) {
    return beta / constants::two_pi_sq;
  }

  //! Convert anisotropic displacement parameters Ucart -> beta.
  /*! This is implemented as Ustar_as_beta(Ucart_as_Ustar(uc, Ucart)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ucart_as_beta(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& Ucart) {
    return Ustar_as_beta(Ucart_as_Ustar(uc, Ucart));
  }
  //! Convert anisotropic displacement parameters beta -> Ucart.
  /*! This is implemented as Ustar_as_Ucart(uc, beta_as_Ustar(beta)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  beta_as_Ucart(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& beta) {
    return Ustar_as_Ucart(uc, beta_as_Ustar(beta));
  }

  //! Convert anisotropic displacement parameters Ucif -> beta.
  /*! This is implemented as Ustar_as_beta(Ucif_as_Ustar(uc, Ucif)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Ucif_as_beta(const uctbx::UnitCell& uc,
               const af::tiny<FloatType, 6>& Ucif) {
    return Ustar_as_beta(Ucif_as_Ustar(uc, Ucif));
  }
  //! Convert anisotropic displacement parameters beta -> Ucif.
  /*! This is implemented as Ustar_as_Ucif(uc, beta_as_Ustar(beta)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  beta_as_Ucif(const uctbx::UnitCell& uc,
               const af::tiny<FloatType, 6>& beta) {
    return Ustar_as_Ucif(uc, beta_as_Ustar(beta));
  }

  //! Convert anisotropic displacement parameters Ucart -> Uiso.
  /*! Uiso is defined as the mean of the diagonal elements of Ucart:<pre>
          Uiso = 1/3 (Ucart11 + Ucart22 + Ucart33)</pre>
   */
  template <typename FloatType>
  inline FloatType
  Ucart_as_Uiso(const af::tiny_plain<FloatType, 6>& Ucart)
  {
    return (Ucart[0] + Ucart[1] + Ucart[2]) / 3.;
  }
  //! Convert Uiso -> anisotropic displacement parameters Ucart.
  /*! The diagonal elements of Ucart are set to the value of Uiso.
      The off-diagonal components Ucart are set to zero.
   */
  template <typename FloatType>
  af::tiny<FloatType, 6>
  Uiso_as_Ucart(const FloatType& Uiso)
  {
    af::tiny<FloatType, 6> result;
    result.fill(0.);
    for(std::size_t i=0;i<3;i++) result[i] = Uiso;
    return result;
  }

  //! Convert Ustar -> Uiso.
  /*! This is implemented as Ucart_as_Uiso(Ustar_as_Ucart(uc, Ustar)).
   */
  template <typename FloatType>
  inline FloatType
  Ustar_as_Uiso(const uctbx::UnitCell& uc,
                const af::tiny<FloatType, 6>& Ustar)
  {
    return Ucart_as_Uiso(Ustar_as_Ucart(uc, Ustar));
  }
  //! Convert Uiso -> Ustar.
  /*! This is implemented as Ucart_as_Ustar(uc, Uiso_as_Ucart(Uiso)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Uiso_as_Ustar(const uctbx::UnitCell& uc,
                const FloatType& Uiso)
  {
    return Ucart_as_Ustar(uc, Uiso_as_Ucart(Uiso));
  }

  //! Convert Ucif -> Uiso.
  /*! This is implemented as Ucart_as_Uiso(Ucif_as_Ucart(uc, Ucif)).
   */
  template <typename FloatType>
  inline FloatType
  Ucif_as_Uiso(const uctbx::UnitCell& uc,
               const af::tiny<FloatType, 6>& Ucif)
  {
    return Ucart_as_Uiso(Ucif_as_Ucart(uc, Ucif));
  }
  //! Convert Uiso -> Ucif.
  /*! This is implemented as Ucart_as_Ucif(uc, Uiso_as_Ucart(Uiso)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  Uiso_as_Ucif(const uctbx::UnitCell& uc,
               const FloatType& Uiso)
  {
    return Ucart_as_Ucif(uc, Uiso_as_Ucart(Uiso));
  }

  //! Convert beta -> Uiso.
  /*! This is implemented as Ucart_as_Uiso(beta_as_Ucart(uc, beta)).
   */
  template <typename FloatType>
  inline FloatType
  beta_as_Uiso(const uctbx::UnitCell& uc,
               const af::tiny<FloatType, 6>& beta)
  {
    return Ucart_as_Uiso(beta_as_Ucart(uc, beta));
  }
  //! Convert Uiso -> beta.
  /*! This is implemented as Ucart_as_beta(uc, Uiso_as_Ucart(Uiso)).
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
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
                        const miller::Index& MIx,
                        double Biso)
  {
    return DebyeWallerFactorBiso(uc.Q(MIx) / 4., Biso);
  }
  //! Isotropic Debye-Waller factor given Miller index and Uiso.
  inline double
  DebyeWallerFactorUiso(const uctbx::UnitCell& uc,
                        const miller::Index& MIx,
                        double Uiso)
  {
    return DebyeWallerFactorBiso(uc, MIx, U_as_B(Uiso));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Ustar.
  template <typename FloatType>
  inline FloatType
  DebyeWallerFactorUstar(const miller::Index& MIx,
                         const af::tiny<FloatType, 6>& Ustar)
  {
    return std::exp(-constants::two_pi_sq * (
        (MIx[0] * MIx[0]) * Ustar[0]
      + (MIx[1] * MIx[1]) * Ustar[1]
      + (MIx[2] * MIx[2]) * Ustar[2]
      + (2 * MIx[0] * MIx[1]) * Ustar[3]
      + (2 * MIx[0] * MIx[2]) * Ustar[4]
      + (2 * MIx[1] * MIx[2]) * Ustar[5]));
  }

  //! Coefficients used in anisotropic Debye-Waller factor calculation.
  /*! Useful for computing partial derivatives w.r.t. Ustar.
   */
  template <typename FloatType>
  inline af::tiny<FloatType, 6>
  DebyeWallerFactorUstarCoefficients(miller::Index const& MIx,
                                     type_holder<FloatType>)
  {
    return -constants::two_pi_sq * af::tiny<FloatType, 6>(
        (MIx[0] * MIx[0]),
        (MIx[1] * MIx[1]),
        (MIx[2] * MIx[2]),
        (2 * MIx[0] * MIx[1]),
        (2 * MIx[0] * MIx[2]),
        (2 * MIx[1] * MIx[2]));
  }

  //! Anisotropic Debye-Waller factor given Miller index and beta.
  template <typename FloatType>
  inline FloatType
  DebyeWallerFactor_beta(const miller::Index& MIx,
                         const af::tiny<FloatType, 6>& beta)
  {
    return DebyeWallerFactorUstar(MIx, beta_as_Ustar(beta));
  }

  //! Anisotropic Debye-Waller factor given Miller index and Ucif.
  template <typename FloatType>
  inline FloatType
  DebyeWallerFactorUcif(const uctbx::UnitCell& uc,
                        const miller::Index& MIx,
                        const af::tiny<FloatType, 6>& Ucif)
  {
    return DebyeWallerFactorUstar(MIx, Ucif_as_Ustar(uc, Ucif));
  }
  //! Anisotropic Debye-Waller factor given Miller index and Ucart.
  template <typename FloatType>
  inline FloatType
  DebyeWallerFactorUcart(const uctbx::UnitCell& uc,
                         const miller::Index& MIx,
                         const af::tiny<FloatType, 6>& Ucart)
  {
    return DebyeWallerFactorUstar(MIx, Ucart_as_Ustar(uc, Ucart));
  }

  namespace detail {

    template <class T>
    inline
    std::complex<T>
    main_root(const std::complex<T>& a, unsigned int m)
    {
      T rho = std::abs(a);
      if (rho == T(0)) return std::complex<T>(0, 0);
      const T one_pi(constants::pi);
      const T two_pi(constants::two_pi);
      T theta0 = std::fmod(std::arg(a), two_pi);
      if (theta0 > one_pi) theta0 -= two_pi;
      if (theta0 <= one_pi) theta0 += two_pi;
      return std::polar(std::pow(rho, T(1./m)), theta0 / T(m));
    }

  } // namespace detail

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
  template <typename FloatType>
  af::tiny<FloatType, 3>
  Eigenvalues(const af::tiny_plain<FloatType, 6>& adp)
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
                  - FloatType(2) * adp[3] * adp[4] * adp[5];
    /* "Reduced form" of the cubic equation according to
       Taschenbuch der Mathematik.
     */
    // reduced form: y^3 + p y + q == 0
    FloatType p = s - r * r / FloatType(3);
    FloatType q =   FloatType(2) * r * r * r / FloatType(27)
                  - r * s / FloatType(3) + t;
    // to circumvent instabilities due to rounding errors the
    // roots are determined as complex numbers.
    std::complex<FloatType> D(p*p*p / FloatType(27) + q*q / FloatType(4));
    std::complex<FloatType> sqrtD = std::sqrt(D);
    FloatType mq2 = -q / FloatType(2);
    std::complex<FloatType> u = detail::main_root(mq2 + sqrtD, 3);
    std::complex<FloatType> v = detail::main_root(mq2 - sqrtD, 3);
    std::complex<FloatType> epsilon1(-0.5, std::sqrt(3.) * 0.5);
    std::complex<FloatType> epsilon2 = std::conj(epsilon1);
    // since the anisotropic displacement tensor is a symmetric matrix,
    // all the imaginary components of the roots must be zero.
    af::tiny<FloatType, 3> result;
    result[0] = (u + v).real();
    result[1] = (epsilon1 * u + epsilon2 * v).real();
    result[2] = (epsilon2 * u + epsilon1 * v).real();
    // convert the solutions y of the reduced form to the
    // solutions x of the normal form.
    for(std::size_t i = 0;i<3;i++) result[i] -= r / FloatType(3);
    return result;
  }

  /*! \brief Test if the anisotropic displacement tensor is
      positive definite, given eigenvalues.
   */
  /*! Test if all eigenvalues are > 0.
      <p>
      See also: CheckPositiveDefinite(), Eigenvalues().
   */
  template <typename FloatType>
  bool
  isPositiveDefinite(const af::tiny<FloatType, 3>& adp_eigenvalues) {
    return adp_eigenvalues[af::min_index(adp_eigenvalues)] > 0.;
  }

  /*! \brief Test if the anisotropic displacement tensor is
      positive definite.
   */
  /*! Test if all eigenvalues are > 0.
      <p>
      See also: CheckPositiveDefinite(), Eigenvalues().
   */
  template <typename FloatType>
  bool
  isPositiveDefinite(const af::tiny<FloatType, 6>& adp) {
    return isPositiveDefinite(Eigenvalues(adp));
  }

  /*! \brief Assert that the anisotropic displacement tensor is
      positive definite, given eigenvalues.
   */
  /*! An exception is thrown if the assertion fails.
      <p>
      See also: isPositiveDefinite(), Eigenvalues().
   */
  template <typename FloatType>
  void
  CheckPositiveDefinite(const af::tiny<FloatType, 3>& adp_eigenvalues) {
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
  template <typename FloatType>
  void
  CheckPositiveDefinite(const af::tiny<FloatType, 6>& adp) {
    CheckPositiveDefinite(Eigenvalues(adp));
  }

  namespace detail {

    template <typename FloatType>
    std::pair<af::tiny<FloatType, 3>, FloatType>
    recursively_multiply(const af::tiny<FloatType, 9>& M,
                         af::tiny<FloatType, 3> V,
                         FloatType tolerance = 1.e-6)
    {
      unsigned int RunAwayCounter = 0;
      for (;;) {
        af::tiny<FloatType, 3> MV;
        MatrixLite::multiply<FloatType>(
          M.begin(), V.begin(), 3, 3, 1, MV.begin());
        FloatType abs_lambda = std::sqrt(af::sum(MV * MV));
        if (abs_lambda == 0.) throw not_positive_definite;
        MV = MV / abs_lambda;
        af::tiny<FloatType, 3> absMV = af::abs(MV);
        std::size_t iMax = af::max_index(absMV);
        FloatType scaled_tolerance = absMV[iMax] * tolerance;
        bool converged = (af::approx_equal_scaled(
          MV, V, scaled_tolerance) == true);
        if (!converged && (af::approx_equal_scaled(
          MV, -V, scaled_tolerance) == true)) {
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
  template <typename FloatType>
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
      Eigensystem(const af::tiny<FloatType, 6>& adp,
                  FloatType tolerance = 1.e-6)
      {
        af::tiny<FloatType, 9> M[2];
        M[0] = MatrixLite::CondensedSymMx33_as_FullSymMx33(adp,
               type_holder<FloatType>());
        FloatType d = MatrixLite::Determinant(M[0]);
        if (d == 0.) {
          throw not_positive_definite;
        }
        M[1] = MatrixLite::CoFactorMxTp(M[0]) / d;
        std::size_t iLarge[2];
        for(std::size_t iM=0;iM<2;iM++) {
          af::tiny<FloatType, 3>
          absDiag = af::abs(MatrixLite::DiagonalElements(M[iM]));
          iLarge[iM] = af::max_index(absDiag);
          if (iM != 0 && iLarge[1] == iLarge[0]) {
            absDiag[iLarge[1]] = -1.;
            iLarge[1] = af::max_index(absDiag);
            cctbx_assert(iLarge[1] != iLarge[0]);
          }
          af::tiny<FloatType, 3> V;
          V.fill(0.);
          V[iLarge[iM]] = 1.;
          std::pair<af::tiny<FloatType, 3>, FloatType>
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
      const af::tiny<FloatType, 3>&
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
      af::tiny<af::tiny<FloatType, 3>, 3> m_vectors;
      af::tiny<FloatType, 3> m_values;
  };

}} // namespace cctbx::adptbx

#endif // CCTBX_ADPTBX_H
