/*! \file
    Toolbox for the handling of anisotropic displacement parameters (ADPs).
 */

#ifndef CCTBX_ADPTBX_H
#define CCTBX_ADPTBX_H

#include <cctbx/uctbx.h>
#include <scitbx/math/eigensystem.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/type_holder.h>
#include <cctbx/error.h>

namespace cctbx {

  //! ADP (anisotropic displacement parameters) Toolbox namespace.
  namespace adptbx {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::sym_mat3;

  inline void
  throw_not_positive_definite()
  {
    throw error("anisotropic displacement tensor is not positive definite.");
  }

  //! Converts isotropic displacement parameter U -> B.
  inline double
  u_as_b(double u_iso)
  {
    return u_iso * scitbx::constants::eight_pi_sq;
  }

  //! Converts isotropic displacement parameter B -> U.
  inline double
  b_as_u(double b_iso)
  {
    return b_iso / scitbx::constants::eight_pi_sq;
  }

  //! Converts anisotropic displacement parameters U -> B.
  template <typename FloatType>
  sym_mat3<FloatType>
  u_as_b(sym_mat3<FloatType> const& u_aniso)
  {
    return
      scitbx::constants::eight_pi_sq * sym_mat3<FloatType>(u_aniso);
  }

  //! Converts anisotropic displacement parameters B -> U.
  template <typename FloatType>
  sym_mat3<FloatType>
  b_as_u(sym_mat3<FloatType> const& b_aniso)
  {
    return FloatType(1. / scitbx::constants::eight_pi_sq) * b_aniso;
  }

  //! Converts anisotropic displacement parameters u_cif -> u_star.
  /*! The transformation matrix used is:<pre>
              (a*  0  0)
          c = ( 0 b*  0)
              ( 0  0 c*)</pre>
      The formula for the transformation is
      u_star = c * u_cif * c.transpose().
      In this particular case the expression simplifies to:<pre>
          u_star_11 = a*^2  u_cif_11
          u_star_22 = b*^2  u_cif_22
          u_star_33 = c*^2  u_cif_33
          u_star_12 = a* b* u_cif_12
          u_star_13 = a* c* u_cif_13
          u_star_23 = b* c* u_cif_23</pre>
   */
  template <typename FloatType>
  sym_mat3<FloatType>
  u_cif_as_u_star(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cif)
  {
    af::double6 const& r_params = unit_cell.reciprocal_parameters();
    return sym_mat3<FloatType>(
      u_cif[0] * (r_params[0] * r_params[0]),
      u_cif[1] * (r_params[1] * r_params[1]),
      u_cif[2] * (r_params[2] * r_params[2]),
      u_cif[3] * (r_params[0] * r_params[1]),
      u_cif[4] * (r_params[0] * r_params[2]),
      u_cif[5] * (r_params[1] * r_params[2]));
  }

  //! Converts anisotropic displacement parameters u_star -> u_cif.
  /*! Inverse of u_cif_as_u_star().
   */
  template <typename FloatType>
  sym_mat3<FloatType>
  u_star_as_u_cif(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_star)
  {
    af::double6 const& r_params = unit_cell.reciprocal_parameters();
    return sym_mat3<FloatType>(
      u_star[0] / (r_params[0] * r_params[0]),
      u_star[1] / (r_params[1] * r_params[1]),
      u_star[2] / (r_params[2] * r_params[2]),
      u_star[3] / (r_params[0] * r_params[1]),
      u_star[4] / (r_params[0] * r_params[2]),
      u_star[5] / (r_params[1] * r_params[2]));
  }

  //! Converts anisotropic displacement parameters u_cart -> u_star.
  /*! The formula for the transformation is
      u_star = c * u_cart * c.transpose(),
      with c = unit_cell.fractionalization_matrix().
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_cart_as_u_star(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cart)
  {
    return u_cart.tensor_transform(unit_cell.fractionalization_matrix());
  }

  //! Converts anisotropic displacement parameters u_star -> u_cart.
  /*! Inverse of u_cart_as_u_star().
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_star_as_u_cart(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_star)
  {
    return u_star.tensor_transform(unit_cell.orthogonalization_matrix());
  }

  //! Converts anisotropic displacement parameters u_cart -> u_cif.
  /*! Implemented without a significant loss of efficiency as
      u_star_as_u_cif(unit_cell, u_cart_as_u_star(unit_cell, u_cart)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_cart_as_u_cif(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cart)
  {
    return u_star_as_u_cif(unit_cell, u_cart_as_u_star(unit_cell, u_cart));
  }

  //! Converts anisotropic displacement parameters u_cif -> u_cart.
  /*! Implemented without a significant loss of efficiency as
      u_star_as_u_cart(unit_cell, u_cif_as_u_star(unit_cell, u_cif)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_cif_as_u_cart(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cif)
  {
    return u_star_as_u_cart(unit_cell, u_cif_as_u_star(unit_cell, u_cif));
  }

  //! Converts anisotropic displacement parameters u_star -> beta.
  /*! The elements of u_star are multiplied by 2pi^2.
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_star_as_beta(sym_mat3<FloatType> const& u_star)
  {
    return FloatType(scitbx::constants::two_pi_sq) * u_star;
  }

  //! Converts anisotropic displacement parameters beta -> u_star.
  /*! The elements of beta are divided by 2pi^2.
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  beta_as_u_star(sym_mat3<FloatType> const& beta)
  {
    return beta / FloatType(scitbx::constants::two_pi_sq);
  }

  //! Converts anisotropic displacement parameters u_cart -> beta.
  /*! Implemented as
      u_star_as_beta(u_cart_as_u_star(unit_cell, u_cart)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_cart_as_beta(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cart)
  {
    return u_star_as_beta(u_cart_as_u_star(unit_cell, u_cart));
  }

  //! Converts anisotropic displacement parameters beta -> u_cart.
  /*! Implemented as
      u_star_as_u_cart(unit_cell, beta_as_u_star(beta)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  beta_as_u_cart(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& beta)
  {
    return u_star_as_u_cart(unit_cell, beta_as_u_star(beta));
  }

  //! Converts anisotropic displacement parameters u_cif -> beta.
  /*! Implemented as
      u_star_as_beta(u_cif_as_u_star(unit_cell, u_cif)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_cif_as_beta(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cif)
  {
    return u_star_as_beta(u_cif_as_u_star(unit_cell, u_cif));
  }

  //! Converts anisotropic displacement parameters beta -> u_cif.
  /*! Implemented as
      u_star_as_u_cif(unit_cell, beta_as_u_star(beta)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  beta_as_u_cif(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& beta)
  {
    return u_star_as_u_cif(unit_cell, beta_as_u_star(beta));
  }

  //! Converts anisotropic displacement parameters u_cart -> u_iso.
  /*! u_iso is defined as the mean of the diagonal elements of u_cart:<pre>
          u_iso = 1/3 (u_cart_11 + u_cart_22 + u_cart_33)</pre>
   */
  template <typename FloatType>
  inline FloatType
  u_cart_as_u_iso(sym_mat3<FloatType> const& u_cart)
  {
    return (u_cart[0] + u_cart[1] + u_cart[2]) / 3.;
  }

  //! Converts u_iso -> anisotropic displacement parameters u_cart.
  /*! The diagonal elements of u_cart are set to the value of u_iso.
      The off-diagonal components u_cart are set to zero.
   */
  template <typename FloatType>
  sym_mat3<FloatType>
  u_iso_as_u_cart(FloatType const& u_iso)
  {
    return sym_mat3<FloatType>(u_iso,u_iso,u_iso,0,0,0);
  }

  //! Converts u_star -> u_iso.
  /*! Implemented as
      u_cart_as_u_iso(u_star_as_u_cart(unit_cell, u_star)).
   */
  template <typename FloatType>
  inline FloatType
  u_star_as_u_iso(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_star)
  {
    return u_cart_as_u_iso(u_star_as_u_cart(unit_cell, u_star));
  }

  //! Converts u_iso -> u_star.
  /*! Implemented as u_cart_as_u_star(unit_cell, u_iso_as_u_cart(u_iso)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_iso_as_u_star(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_iso)
  {
    return u_cart_as_u_star(unit_cell, u_iso_as_u_cart(u_iso));
  }

  //! Converts u_cif -> u_iso.
  /*! Implemented as
      u_cart_as_u_iso(u_cif_as_u_cart(unit_cell, u_cif)).
   */
  template <typename FloatType>
  inline FloatType
  u_cif_as_u_iso(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& u_cif)
  {
    return u_cart_as_u_iso(u_cif_as_u_cart(unit_cell, u_cif));
  }

  //! Converts u_iso -> u_cif.
  /*! Implemented as
      u_cart_as_u_cif(unit_cell, u_iso_as_u_cart(u_iso)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_iso_as_u_cif(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_iso)
  {
    return u_cart_as_u_cif(unit_cell, u_iso_as_u_cart(u_iso));
  }

  //! Converts beta -> u_iso.
  /*! Implemented as
      u_cart_as_u_iso(beta_as_u_cart(unit_cell, beta)).
   */
  template <typename FloatType>
  inline FloatType
  beta_as_u_iso(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& beta)
  {
    return u_cart_as_u_iso(beta_as_u_cart(unit_cell, beta));
  }

  //! Converts u_iso -> beta.
  /*! Implemented as
      u_cart_as_beta(unit_cell, u_iso_as_u_cart(u_iso)).
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  u_iso_as_beta(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_iso)
  {
    return u_cart_as_beta(unit_cell, u_iso_as_u_cart(u_iso));
  }

  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and b_iso.
  inline double
  debye_waller_factor_b_iso(
    double stol_sq,
    double b_iso)
  {
    return std::exp(-b_iso * stol_sq);
  }

  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and u_iso.
  inline double
  debye_waller_factor_u_iso(
    double stol_sq,
    double u_iso)
  {
    return debye_waller_factor_b_iso(stol_sq, u_as_b(u_iso));
  }

  //! Isotropic Debye-Waller factor given a Miller index and b_iso.
  inline double
  debye_waller_factor_b_iso(
    uctbx::unit_cell const& unit_cell,
    miller::index<> const& h,
    double b_iso)
  {
    return debye_waller_factor_b_iso(unit_cell.d_star_sq(h) / 4., b_iso);
  }

  //! Isotropic Debye-Waller factor given a Miller index and u_iso.
  inline double
  debye_waller_factor_u_iso(
    uctbx::unit_cell const& unit_cell,
    miller::index<> const& h,
    double u_iso)
  {
    return debye_waller_factor_b_iso(unit_cell, h, u_as_b(u_iso));
  }

  //! Anisotropic Debye-Waller factor given a Miller index and u_star.
  template <typename FloatType>
  inline FloatType
  debye_waller_factor_u_star(
    miller::index<> const& h,
    sym_mat3<FloatType> const& u_star)
  {
    return std::exp(-scitbx::constants::two_pi_sq * (
        (h[0] * h[0]) * u_star[0]
      + (h[1] * h[1]) * u_star[1]
      + (h[2] * h[2]) * u_star[2]
      + (2 * h[0] * h[1]) * u_star[3]
      + (2 * h[0] * h[2]) * u_star[4]
      + (2 * h[1] * h[2]) * u_star[5]));
  }

  //! Coefficients used in anisotropic Debye-Waller factor calculation.
  /*! Useful for computing partial gradients w.r.t. u_star.
      <p>
      Not available in Python.
   */
  template <typename NumType>
  inline sym_mat3<NumType>
  debye_waller_factor_u_star_coefficients(
    miller::index<> const& h,
    scitbx::type_holder<NumType>)
  {
    return sym_mat3<NumType>(
        (h[0] * h[0]),
        (h[1] * h[1]),
        (h[2] * h[2]),
        (2 * h[0] * h[1]),
        (2 * h[0] * h[2]),
        (2 * h[1] * h[2]));
  }

  //! Anisotropic Debye-Waller factor given a Miller index and beta.
  template <typename FloatType>
  inline FloatType
  debye_waller_factor_beta(
    miller::index<> const& h,
    sym_mat3<FloatType> const& beta)
  {
    return debye_waller_factor_u_star(h, beta_as_u_star(beta));
  }

  //! Anisotropic Debye-Waller factor given a Miller index and u_cif.
  template <typename FloatType>
  inline FloatType
  debye_waller_factor_u_cif(
    uctbx::unit_cell const& unit_cell,
    miller::index<> const& h,
    sym_mat3<FloatType> const& u_cif)
  {
    return debye_waller_factor_u_star(h, u_cif_as_u_star(unit_cell, u_cif));
  }

  //! Anisotropic Debye-Waller factor given a Miller index and u_cart.
  template <typename FloatType>
  inline FloatType
  debye_waller_factor_u_cart(
    uctbx::unit_cell const& unit_cell,
    miller::index<> const& h,
    sym_mat3<FloatType> const& u_cart)
  {
    return debye_waller_factor_u_star(h, u_cart_as_u_star(unit_cell, u_cart));
  }

  namespace detail {

    /* Mathematica script used to determine the transformation law:
         SetOptions["stdout", PageWidth -> 50]
         us={{s00,s01,s02},{s01,s11,s12},{s02,s12,s22}}
         hs={hs0,hs1,hs2}
         fs=Exp[cb hs.us.hs]
         gs={{D[fs,s00],D[fs,s01],D[fs,s02]},
             {D[fs,s01],D[fs,s11],D[fs,s12]},
             {D[fs,s02],D[fs,s12],D[fs,s22]}}/fs
         uc={{c00,c01,c02},{c01,c11,c12},{c02,c12,c22}}
         hc={hc0,hc1,hc2}
         fc=Exp[cb hc.uc.hc]
         gc={{D[fc,c00],D[fc,c01],D[fc,c02]},
             {D[fc,c01],D[fc,c11],D[fc,c12]},
             {D[fc,c02],D[fc,c12],D[fc,c22]}}/fc
         a={{a00,a01,a02},
            {a10,a11,a12},
            {a20,a21,a22}}
         hc=Transpose[a].hs
         hc0=hc[[1]]
         hc1=hc[[2]]
         hc2=hc[[3]]
         FullSimplify[gs]
         gcs=Expand[Expand[FullSimplify[gc]]
           /. cb hs0^2 -> g00
           /. cb hs1^2 -> g11
           /. cb hs2^2 -> g22
           /. cb hs0 hs1 -> g01/2
           /. cb hs0 hs2 -> g02/2
           /. cb hs1 hs2 -> g12/2]
         InputForm[gcs[[1,1]]]
         InputForm[gcs[[2,2]]]
         InputForm[gcs[[3,3]]]
         InputForm[gcs[[1,2]]]
         InputForm[gcs[[1,3]]]
         InputForm[gcs[[2,3]]]
     */
    template <typename FloatType>
    inline sym_mat3<FloatType>
    grad_u_transform(
      mat3<FloatType> const& a,
      sym_mat3<FloatType> const& g)
    {
      return sym_mat3<FloatType>(
          a[0]*a[0]*g[0] + a[0]*a[3]*g[3] + a[0]*a[6]*g[4]
        + a[3]*a[3]*g[1] + a[3]*a[6]*g[5] + a[6]*a[6]*g[2],
          a[1]*a[1]*g[0] + a[1]*a[4]*g[3] + a[1]*a[7]*g[4]
        + a[4]*a[4]*g[1] + a[4]*a[7]*g[5] + a[7]*a[7]*g[2],
          a[2]*a[2]*g[0] + a[2]*a[5]*g[3] + a[2]*a[8]*g[4]
        + a[5]*a[5]*g[1] + a[5]*a[8]*g[5] + a[8]*a[8]*g[2],
          2*a[0]*a[1]*g[0] + a[1]*a[3]*g[3] +   a[0]*a[4]*g[3]
        +   a[1]*a[6]*g[4] + a[0]*a[7]*g[4] + 2*a[3]*a[4]*g[1]
        +   a[4]*a[6]*g[5] + a[3]*a[7]*g[5] + 2*a[6]*a[7]*g[2],
          2*a[0]*a[2]*g[0] + a[2]*a[3]*g[3] +   a[0]*a[5]*g[3]
        +   a[2]*a[6]*g[4] + a[0]*a[8]*g[4] + 2*a[3]*a[5]*g[1]
        +   a[5]*a[6]*g[5] + a[3]*a[8]*g[5] + 2*a[6]*a[8]*g[2],
          2*a[1]*a[2]*g[0] + a[2]*a[4]*g[3] +   a[1]*a[5]*g[3]
        +   a[2]*a[7]*g[4] + a[1]*a[8]*g[4] + 2*a[4]*a[5]*g[1]
        +   a[5]*a[7]*g[5] + a[4]*a[8]*g[5] + 2*a[7]*a[8]*g[2]);
    }

  } // namespace detail

  //! Transformation of gradients w.r.t. u_star to gradients w.r.t. u_cart.
  /*! Scalar version.
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  grad_u_star_as_u_cart(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& grad_u_star)
  {
    return detail::grad_u_transform(
      unit_cell.fractionalization_matrix(), grad_u_star);
  }

  //! Transformation of gradients w.r.t. u_star to gradients w.r.t. u_cart.
  /*! Vector version.
   */
  template <typename FloatType>
  af::shared<sym_mat3<FloatType> >
  grad_u_star_as_u_cart(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<sym_mat3<FloatType> > const& grad_u_star)
  {
    af::shared<sym_mat3<FloatType> > result((af::reserve(grad_u_star.size())));
    for(std::size_t i=0;i<grad_u_star.size();i++) {
      result.push_back(grad_u_star_as_u_cart(unit_cell, grad_u_star[i]));
    }
    return result;
  }

  //! Transformation of gradients w.r.t. u_cart to gradients w.r.t. u_star.
  /*! Scalar version.
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  grad_u_cart_as_u_star(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& grad_u_cart)
  {
    return detail::grad_u_transform(
      unit_cell.orthogonalization_matrix(), grad_u_cart);
  }

  //! Transformation of gradients w.r.t. u_cart to gradients w.r.t. u_star.
  /*! Vector version.
   */
  template <typename FloatType>
  af::shared<sym_mat3<FloatType> >
  grad_u_cart_as_u_star(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<sym_mat3<FloatType> > const& grad_u_cart)
  {
    af::shared<sym_mat3<FloatType> > result((af::reserve(grad_u_cart.size())));
    for(std::size_t i=0;i<grad_u_cart.size();i++) {
      result.push_back(grad_u_cart_as_u_star(unit_cell, grad_u_cart[i]));
    }
    return result;
  }

  namespace detail {

    template <class T>
    inline
    std::complex<T>
    root(std::complex<T> const& c, unsigned int m, unsigned int k=0)
    {
      T rho = std::abs(c);
      if (rho == T(0)) return std::complex<T>(0, 0);
      return std::polar(
        std::pow(rho, T(1./m)),
        (std::arg(c) + T(k * scitbx::constants::two_pi)) / T(m));
    }

    // To ensure numerical stability we have to try all 3 roots
    // of one of the complex numbers. The goal is to find the root
    // which minimizes the sum of the imaginary parts of u and v.
    template <typename FloatType>
    struct third_roots
    {
      third_roots(
        std::complex<FloatType> const& c1,
        std::complex<FloatType> const& c2)
      {
        u = detail::root(c1, 3);
        v = detail::root(c2, 3);
        FloatType diff = scitbx::fn::absolute(u.imag() + v.imag());
        for(unsigned int k=1;k<3;k++) {
          std::complex<FloatType> trial_v = detail::root(c2, 3, k);
          FloatType trial_diff=scitbx::fn::absolute(u.imag() + trial_v.imag());
          if (diff > trial_diff) {
            diff = trial_diff;
            v = trial_v;
          }
        }
      }

      std::complex<FloatType> u;
      std::complex<FloatType> v;
    };

  } // namespace detail

  //! Determines the eigenvalues of the anisotropic displacement tensor.
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
      See also: eigensystem().
   */
  template <typename FloatType>
  vec3<FloatType>
  eigenvalues(sym_mat3<FloatType> const& adp)
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
    FloatType d(p*p*p / FloatType(27) + q*q / FloatType(4));
    std::complex<FloatType> sqrt_d;
    if (d < FloatType(0)) {
      sqrt_d = std::complex<FloatType>(0, std::sqrt(-d));
    }
    else {
      sqrt_d = std::complex<FloatType>(std::sqrt(d), 0);
    }
    FloatType mq2 = -q / FloatType(2);
    detail::third_roots<FloatType> roots(mq2 + sqrt_d, mq2 - sqrt_d);
    std::complex<FloatType> epsilon1(-0.5, std::sqrt(3.) * 0.5);
    std::complex<FloatType> epsilon2 = std::conj(epsilon1);
    // since the anisotropic displacement tensor is a symmetric matrix,
    // all the imaginary components of the roots must be zero.
    vec3<FloatType> result;
    result[0] = (roots.u + roots.v).real();
    result[1] = (epsilon1 * roots.u + epsilon2 * roots.v).real();
    result[2] = (epsilon2 * roots.u + epsilon1 * roots.v).real();
    // convert the solutions y of the reduced form to the
    // solutions x of the normal form.
    for(std::size_t i = 0;i<3;i++) result[i] -= r / FloatType(3);
    return result;
  }

  /*! \brief Tests if the anisotropic displacement tensor is
      positive definite, given adp_eigenvalues.
   */
  /*! Tests if all adp_eigenvalues are > 0.
      <p>
      See also: eigenvalues().
   */
  template <typename FloatType>
  bool
  is_positive_definite(
    vec3<FloatType> const& adp_eigenvalues)
  {
    return scitbx::af::min(adp_eigenvalues.const_ref()) > 0;
  }

  /*! \brief Tests if the anisotropic displacement tensor is
      positive definite, given adp_eigenvalues.
   */
  /*! Tests if all adp_eigenvalues are >= -tolerance.
      <p>
      See also: eigenvalues().
   */
  template <typename FloatType>
  bool
  is_positive_definite(
    vec3<FloatType> const& adp_eigenvalues,
    FloatType const& tolerance)
  {
    return scitbx::af::min(adp_eigenvalues.const_ref()) >= -tolerance;
  }

  /*! \brief Tests if the anisotropic displacement tensor is
      positive definite.
   */
  /*! Tests if all eigenvalues(adp) are > 0.
   */
  template <typename FloatType>
  bool
  is_positive_definite(
    sym_mat3<FloatType> const& adp)
  {
    return is_positive_definite(eigenvalues(adp));
  }

  /*! \brief Tests if the anisotropic displacement tensor is
      positive definite.
   */
  /*! Tests if all eigenvalues(adp) are >= -tolerance.
   */
  template <typename FloatType>
  bool
  is_positive_definite(
    sym_mat3<FloatType> const& adp,
    FloatType const& tolerance)
  {
    return is_positive_definite(eigenvalues(adp), tolerance);
  }

  namespace detail {

    template <typename FloatType>
    struct recursively_multiply
    {
      vec3<FloatType> v;
      FloatType abs_lambda;

      recursively_multiply(mat3<FloatType> const& m,
                           vec3<FloatType> const& v_start,
                           FloatType tolerance=1.e-6)
      :
        v(v_start)
      {
        unsigned int run_away_counter = 0;
        for (;;) {
          vec3<FloatType> mv = m * v;
          abs_lambda = mv.length();
          if (abs_lambda == 0.) throw_not_positive_definite();
          mv /= abs_lambda;
          vec3<FloatType> abs_mv = af::fabs(mv.as_tiny());
          std::size_t i_max = af::max_index(abs_mv);
          FloatType scaled_tolerance = abs_mv[i_max] * tolerance;
          bool converged = mv.const_ref().all_approx_equal(
                            v.const_ref(), scaled_tolerance);
          if (!converged && mv.const_ref().all_approx_equal(
                            (-v).const_ref(), scaled_tolerance)) {
            throw_not_positive_definite(); // lambda < 0
          }
          v = mv;
          if (converged) return;
          run_away_counter++;
          CCTBX_ASSERT(run_away_counter < 10000000);
        }
      }
    };

  } // namespace detail

  //! Group of associated eigenvectors and eigenvalues.
  template <typename FloatType>
  class eigensystem
  {
    public:
      //! Default constructor. Some data members are not initialized!
      eigensystem() {}

      /*! \brief Determines the eigenvectors and eigenvalues of the
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
          See also: eigenvalues().
       */
      eigensystem(sym_mat3<FloatType> const& adp,
                  FloatType tolerance=1.e-6);

      //! The i'th eigenvector.
      /*! An exception is thrown if i >= 3.
       */
      vec3<FloatType> const&
      vectors(std::size_t i) const
      {
        if (i >= vectors_.size()) throw error_index();
        return vectors_[i];
      }

      //! The eigenvalues.
      vec3<FloatType> const&
      values() const { return values_; }

    private:
      af::tiny<vec3<FloatType>, 3> vectors_;
      vec3<FloatType> values_;
  };

  template <typename FloatType>
  eigensystem<FloatType>
  ::eigensystem(sym_mat3<FloatType> const& adp,
                FloatType tolerance)
  {
    mat3<FloatType> m[2];
    m[0] = mat3<FloatType>(adp);
    FloatType d = m[0].determinant();
    if (d == 0.) throw_not_positive_definite();
    m[1] = m[0].co_factor_matrix_transposed() / d;
    std::size_t i_large[2];
    for(std::size_t i_m=0;i_m<2;i_m++) {
      vec3<FloatType> abs_diag = af::fabs(m[i_m].diagonal().as_tiny());
      i_large[i_m] = af::max_index(abs_diag);
      if (i_m != 0 && i_large[1] == i_large[0]) {
        abs_diag[i_large[1]] = -1.;
        i_large[1] = af::max_index(abs_diag);
        CCTBX_ASSERT(i_large[1] != i_large[0]);
      }
      vec3<FloatType> v(0,0,0);
      v[i_large[i_m]] = 1.;
      detail::recursively_multiply<FloatType> v_lambda(m[i_m], v);
      vectors_[i_m*2] = v_lambda.v;
      values_[i_m*2] = v_lambda.abs_lambda;
    }
    vectors_[1] = vectors_[0].cross(vectors_[2]);
    CCTBX_ASSERT(vectors_[1].length() != 0);
    values_[2] = 1. / values_[2];
    values_[1] = (adp[0] + adp[1] + adp[2]) - (values_[0] + values_[2]);
  }

  //! Modifies u_cart such that all eigenvalues are >= 0.
  template <typename FloatType>
  sym_mat3<FloatType>
  eigenvalue_filtering(sym_mat3<FloatType> const& u_cart)
  {
    scitbx::math::eigensystem::real_symmetric<FloatType> es(u_cart);
    scitbx::vec3<FloatType> es_val(es.values().begin());
    for(std::size_t i=0;i<3;i++) if (es_val[i] < 0) es_val[i] = 0;
    scitbx::mat3<FloatType> es_vec(es.vectors().begin());
    scitbx::mat3<FloatType> es_vec_inv = es_vec.inverse();
    return sym_mat3<FloatType>(es_val).tensor_transform(es_vec_inv);
  }

  //! Tensor transformation: c * u * c.transpose().
  /*! For use in Python only.
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  c_u_c_transpose(mat3<FloatType> const& c, sym_mat3<FloatType> const& u)
  {
    return u.tensor_transform(c);
  }

}} // namespace cctbx::adptbx

#endif // CCTBX_ADPTBX_H
