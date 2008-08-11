#ifndef CCTBX_ELTBX_XRAY_SCATTERING_GAUSSIAN_H
#define CCTBX_ELTBX_XRAY_SCATTERING_GAUSSIAN_H

#include <cctbx/eltbx/xray_scattering/form_factor.h>
#include <scitbx/math/gaussian/sum.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

  class gaussian
    : public scitbx::math::gaussian::sum<double>,
      public isotropic_form_factor_mixin<gaussian>
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
  };

}}} // cctbx::eltbx::xray_scattering

#endif // CCTBX_ELTBX_XRAY_SCATTERING_GAUSSIAN_H
