/*! \file
    Toolbox for the handling of anisotropic displacement parameters (ADPs).
 */

#ifndef CCTBX_ADPTBX_H
#define CCTBX_ADPTBX_H

#include <cctbx/uctbx.h>
#include <scitbx/matrix/eigensystem.h>
#include <scitbx/matrix/tensor_rank_2.h>
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

  //! Factorization into isotropic and remaining anisotropic contributions.
  template <typename FloatType=double>
  struct factor_u_cart_u_iso
  {
    factor_u_cart_u_iso() {}

    factor_u_cart_u_iso(sym_mat3<FloatType> const& u_cart)
    :
      u_iso(u_cart_as_u_iso(u_cart)),
      u_cart_minus_u_iso(u_cart)
    {
      for(unsigned i=0;i<3;i++) {
        u_cart_minus_u_iso[i] -= u_iso;
      }
    }

    FloatType u_iso;
    sym_mat3<FloatType> u_cart_minus_u_iso;
  };

  //! Factorization into isotropic and remaining anisotropic contributions.
  template <typename FloatType=double>
  struct factor_u_star_u_iso
  {
    factor_u_star_u_iso() {}

    factor_u_star_u_iso(
      uctbx::unit_cell const& unit_cell,
      sym_mat3<FloatType> const& u_star)
    {
      factor_u_cart_u_iso<FloatType>
        f_cart(u_star_as_u_cart(unit_cell, u_star));
      u_iso = f_cart.u_iso;
      u_star_minus_u_iso = u_cart_as_u_star(
        unit_cell, f_cart.u_cart_minus_u_iso);
    }

    FloatType u_iso;
    sym_mat3<FloatType> u_star_minus_u_iso;
  };

  //! Factorization into isotropic and remaining anisotropic contributions.
  template <typename FloatType=double>
  struct factor_u_cif_u_iso
  {
    factor_u_cif_u_iso() {}

    factor_u_cif_u_iso(
      uctbx::unit_cell const& unit_cell,
      sym_mat3<FloatType> const& u_cif)
    {
      factor_u_cart_u_iso<FloatType>
        f_cart(u_cif_as_u_cart(unit_cell, u_cif));
      u_iso = f_cart.u_iso;
      u_cif_minus_u_iso = u_cart_as_u_cif(
        unit_cell, f_cart.u_cart_minus_u_iso);
    }

    FloatType u_iso;
    sym_mat3<FloatType> u_cif_minus_u_iso;
  };

  //! Factorization into isotropic and remaining anisotropic contributions.
  template <typename FloatType=double>
  struct factor_beta_u_iso
  {
    factor_beta_u_iso() {}

    factor_beta_u_iso(
      uctbx::unit_cell const& unit_cell,
      sym_mat3<FloatType> const& beta)
    {
      factor_u_cart_u_iso<FloatType>
        f_cart(beta_as_u_cart(unit_cell, beta));
      u_iso = f_cart.u_iso;
      beta_minus_u_iso = u_cart_as_beta(
        unit_cell, f_cart.u_cart_minus_u_iso);
    }

    FloatType u_iso;
    sym_mat3<FloatType> beta_minus_u_iso;
  };

  //! std::exp with upper limit for argument value.
  inline double
  debye_waller_factor_exp(const char* u_type, double arg, double max_arg=50)
  {
    if (arg > max_arg) {
      char buf[256];
      std::sprintf(buf,
        "cctbx::adptbx::debye_waller_factor_exp:"
        " max_arg exceeded (%s):"
        " arg = %.6g"
        " max_arg = %.6g",
        u_type, arg, max_arg);
      throw std::runtime_error(buf);
    }
    return std::exp(arg);
  }

  //! Isotropic Debye-Waller factor given (sin(theta)/lambda)^2 and b_iso.
  inline double
  debye_waller_factor_b_iso(
    double stol_sq,
    double b_iso)
  {
    return debye_waller_factor_exp("isotropic", -b_iso * stol_sq);
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
    return debye_waller_factor_exp(
      "anisotropic",
      -scitbx::constants::two_pi_sq * (
          (h[0] * h[0]) * u_star[0]
        + (h[1] * h[1]) * u_star[1]
        + (h[2] * h[2]) * u_star[2]
        + (2 * h[0] * h[1]) * u_star[3]
        + (2 * h[0] * h[2]) * u_star[4]
        + (2 * h[1] * h[2]) * u_star[5]));
  }

  //! Coefficients for gradients of Debye-Waller factor w.r.t. u_star.
  /*! Formula for the gradients:
        -2*pi**2 * debye_waller_factor_u_star(h, u_star) * result

      It should be noted that passing a 2nd argument, e.g.
      \code
        debye_waller_factor_u_star_gradient_coefficients(
          h, scitbx::type_holder<double>());
      \endcode
      to select the version returning sym_mat3<double>, is deprecated.
      All the platform supported by the cctbx have C++ compiler which accepts
      the simpler syntax
      \code
        debye_waller_factor_u_star_gradient_coefficients<double>(h)
      \endcode
      The last one which failed to support such function template
      specialisation, VC++ 7.1, has been dropped in early 2009.
   */
  template <typename NumType>
  inline sym_mat3<NumType>
  debye_waller_factor_u_star_gradient_coefficients(
    miller::index<> const& h,
    scitbx::type_holder<NumType> holder=scitbx::type_holder<NumType>())
  {
    return sym_mat3<NumType>(
        (h[0] * h[0]),
        (h[1] * h[1]),
        (h[2] * h[2]),
        (2 * h[0] * h[1]),
        (2 * h[0] * h[2]),
        (2 * h[1] * h[2]));
  }

  //! Coefficients for curvatures of Debye-Waller factor w.r.t. u_star.
  /*! Formula for the curvatures:
        (-2*pi**2)**2 * debye_waller_factor_u_star(h, u_star) * result

      The returned array contains the 6*(6+1)/2 elements of the
      upper diagonal of the (6 x 6) matrix of curvatures.
   */
  template <typename NumType>
  af::shared<NumType>
  debye_waller_factor_u_star_curvature_coefficients(
    miller::index<> const& h,
    scitbx::type_holder<NumType>)
  {
    /* Python script for generating the code below:

         sym_mat3_indices = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
         i = 0
         for ij1 in xrange(6):
           i1,j1 = sym_mat3_indices[ij1]
           for ij2 in xrange(ij1,6):
             i2,j2 = sym_mat3_indices[ij2]
             print "result[%d] = h%dh%d * h%dh%d;" % (i,i1,j1,i2,j2)
             i += 1
     */
    af::shared<NumType> result(6*(6+1)/2, 0);
    NumType h0h0 = h[0] * h[0];
    NumType h1h1 = h[1] * h[1];
    NumType h2h2 = h[2] * h[2];
    NumType h0h1 = 2 * h[0] * h[1];
    NumType h0h2 = 2 * h[0] * h[2];
    NumType h1h2 = 2 * h[1] * h[2];
    result[0] = h0h0 * h0h0;
    result[1] = h0h0 * h1h1;
    result[2] = h0h0 * h2h2;
    result[3] = h0h0 * h0h1;
    result[4] = h0h0 * h0h2;
    result[5] = h0h0 * h1h2;
    result[6] = h1h1 * h1h1;
    result[7] = h1h1 * h2h2;
    result[8] = h1h1 * h0h1;
    result[9] = h1h1 * h0h2;
    result[10] = h1h1 * h1h2;
    result[11] = h2h2 * h2h2;
    result[12] = h2h2 * h0h1;
    result[13] = h2h2 * h0h2;
    result[14] = h2h2 * h1h2;
    result[15] = h0h1 * h0h1;
    result[16] = h0h1 * h0h2;
    result[17] = h0h1 * h1h2;
    result[18] = h0h2 * h0h2;
    result[19] = h0h2 * h1h2;
    result[20] = h1h2 * h1h2;
    return result;
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

  //! Transformation of gradients w.r.t. u_star to gradients w.r.t. u_cart.
  /*! Scalar version.
   */
  template <typename FloatType>
  inline sym_mat3<FloatType>
  grad_u_star_as_u_cart(
    uctbx::unit_cell const& unit_cell,
    sym_mat3<FloatType> const& grad_u_star)
  {
    return scitbx::matrix::tensor_rank_2::gradient_transform(
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
    return scitbx::matrix::tensor_rank_2::gradient_transform(
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

          This class is implemented as a thin wrapper around:
            scitbx::matrix::eigensystem::real_symmetric<>
       */
      eigensystem(sym_mat3<FloatType> const& adp)
      {
        scitbx::matrix::eigensystem::real_symmetric<FloatType> es(adp);
        for(std::size_t i=0;i<3;i++) {
          vectors_[i] = vec3<FloatType>(&es.vectors()[i*3]);
        }
        values_ = vec3<FloatType>(es.values().begin());
      }

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

  //! Determines the eigenvalues of the anisotropic displacement tensor.
  /*! Equivalent to eigensystem<>().values().
   */
  template <typename FloatType>
  vec3<FloatType>
  eigenvalues(sym_mat3<FloatType> const& adp)
  {
    return eigensystem<FloatType>(adp).values();
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

  /*! \brief True if adp is positive definite, else False.
   */
  /*! Tests if all eigenvalues(adp) are >= -tolerance.
   */
  template <typename FloatType>
  af::shared<bool>
  is_positive_definite(
    af::const_ref<sym_mat3<FloatType> > const& adp,
    FloatType const& tolerance)
  {
    af::shared<bool> result((af::reserve(adp.size())));
    for(std::size_t i=0;i<adp.size();i++) {
        result.push_back(is_positive_definite(eigenvalues(adp[i]), tolerance));
    }
    return result;
  }

  //! Isotropize u_cart: modify u_cart such that it meets target anisotropy.
  template <typename FloatType>
  sym_mat3<FloatType>
  isotropize(
    sym_mat3<FloatType> const& u_cart,
    FloatType const& anisotropy_min=0.25)
  {
    scitbx::matrix::eigensystem::real_symmetric<FloatType> es(u_cart);
    scitbx::vec3<FloatType> es_val(es.values().begin());
    FloatType u_min = es_val[0];
    FloatType u_max = u_min;
    for(std::size_t i=0;i<3;i++) {
      if(es_val[i] < u_min) u_min = es_val[i];
      if(es_val[i] > u_max) u_max = es_val[i];
    }
    FloatType anisotropy = anisotropy_min;
    if(u_max != 0) anisotropy = u_min/u_max;
    FloatType corr = 0;
    if(anisotropy < anisotropy_min) {
      corr = (anisotropy_min * u_max - u_min)/(anisotropy_min + 1);
      for(std::size_t i=0;i<3;i++) {
        if(es_val[i] == u_min) es_val[i] = u_min + corr;
        if(es_val[i] == u_max) es_val[i] = u_max - corr;
      }
      scitbx::mat3<FloatType> es_vec(es.vectors().begin());
      scitbx::mat3<FloatType> es_vec_inv = es_vec.inverse();
      return sym_mat3<FloatType>(es_val).tensor_transform(es_vec_inv);
    }
    else return u_cart;
  }

  //! Modifies u_cart such that all eigenvalues are >= u_min and <= u_max.
  /*! u_max is used only if it is greater than zero.
   */
  template <typename FloatType>
  sym_mat3<FloatType>
  eigenvalue_filtering(
    sym_mat3<FloatType> const& u_cart,
    FloatType const& u_min=0,
    FloatType const& u_max=0)
  {
    scitbx::matrix::eigensystem::real_symmetric<FloatType> es(u_cart);
    scitbx::vec3<FloatType> es_val(es.values().begin());
    for(std::size_t i=0;i<3;i++) if (es_val[i] < u_min) es_val[i] = u_min;
    if (u_max > 0) {
      for(std::size_t i=0;i<3;i++) if (es_val[i] > u_max) es_val[i] = u_max;
    }
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
