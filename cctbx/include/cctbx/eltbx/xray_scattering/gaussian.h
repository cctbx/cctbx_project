#ifndef CCTBX_ELTBX_XRAY_SCATTERING_GAUSSIAN_H
#define CCTBX_ELTBX_XRAY_SCATTERING_GAUSSIAN_H

#include <scitbx/math/gaussian/sum.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

  class gaussian : public scitbx::math::gaussian::sum<double>
  {
    public:
      typedef scitbx::math::gaussian::sum<double> base_t;

      //! Default constructor. Some data members are not initialized!
      gaussian() {}

      //! Initialization given an instance of the base type.
      gaussian(base_t const& gaussian_sum)
      :
        base_t(gaussian_sum)
      {}

      //! Initialization of the constant.
      explicit
      gaussian(double c, bool use_c=true)
      :
        base_t(c, use_c)
      {}

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true.
       */
      gaussian(
        af::small<double, base_t::max_n_terms> const& a,
        af::small<double, base_t::max_n_terms> const& b,
        double c=0,
        bool use_c=false)
      :
        base_t(a, b, c, use_c)
      {}

      //! Initialization of the terms and optionally the constant.
      /*! If c is different from zero use_c will automatically be
          set to true.
       */
      gaussian(
        af::const_ref<double> const& ab,
        double c=0,
        bool use_c=false)
      :
        base_t(ab, c, use_c)
      {}

      /*! \brief Sum of Gaussian terms at the point stol
          (sin-theta-over-lambda), given stol^2.
       */
      /*! See also: at_stol(), at_d_star(), at_d_star_sq(),
                    uctbx::unit_cell::stol()
       */
      double
      at_stol_sq(double stol_sq) const
      {
        return at_x_sq(stol_sq);
      }

      /*! \brief Sum of Gaussian terms at the point stol
          (sin-theta-over-lambda).
       */
      /*! See also: at_stol_sq(), at_d_star(), at_d_star_sq(),
       */
      double
      at_stol(double stol) const
      {
        return at_x_sq(stol * stol);
      }

      /*! \brief Sum of Gaussian terms at the point d_star
          (1/d), given d_star^2.
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star(),
                    uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star_sq(double d_star_sq) const
      {
        return at_x_sq(d_star_sq / 4);
      }

      /*! \brief Sum of Gaussian terms at the point d_star
          (1/d).
       */
      /*! See also: at_stol_sq(), at_stol(), at_d_star_sq(),
                    uctbx::unit_cell::d_star_sq()
       */
      double
      at_d_star(double d_star) const
      {
        return at_x_sq(d_star * d_star / 4);
      }
  };

}}} // cctbx::eltbx::xray_scattering

#endif // CCTBX_ELTBX_XRAY_SCATTERING_GAUSSIAN_H
