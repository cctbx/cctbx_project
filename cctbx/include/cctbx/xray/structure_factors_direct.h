#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/sgtbx/miller_ops.h>
#include <cctbx/math/cos_sin_table.h>

namespace cctbx { namespace xray { namespace structure_factors {

  template <typename CosSinType, typename ScattererType>
  struct direct_one_h_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    direct_one_h_one_scatterer(
      CosSinType const& cos_sin,
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      float_type d_star_sq,
      ScattererType const& scatterer)
    :
      f_calc(0,0)
    {
      typedef float_type f_t;
      typedef std::complex<f_t> c_t;
      for(std::size_t i_smx=0;i_smx<space_group.n_smx();i_smx++) {
        sgtbx::rt_mx const& s = space_group.smx(i_smx);
        miller::index<> hr = h * s.r();
        f_t hrx = hr * scatterer.site;
        f_t ht = f_t(h * s.t()) / space_group.t_den();
        c_t term = cos_sin.get(hrx + ht);
        if (scatterer.anisotropic_flag) {
          f_t dw = adptbx::debye_waller_factor_u_star(hr, scatterer.u_star);
          term *= dw;
        }
        f_calc += term;
      }
      if (!scatterer.anisotropic_flag && scatterer.u_iso != 0) {
        f_t dw=adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        f_calc *= dw;
      }
    }

    std::complex<float_type> f_calc;
  };

  template <typename CosSinType, typename ScattererType>
  struct direct_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    direct_one_scatterer(
      CosSinType const& cos_sin,
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      af::const_ref<miller::index<> > const& miller_indices,
      ScattererType const& scatterer,
      af::ref<std::complex<float_type> > const& f_calc)
    {
      typedef float_type f_t;
      typedef std::complex<float_type> c_t;
      f_t w = scatterer.weight();
      f_t fp_w = scatterer.fp;
      f_t fdp_w = scatterer.fdp;
      bool have_fdp = fdp_w != 0;
      bool caasf_is_const = std::strcmp(scatterer.caasf.label(), "const") == 0;
      if (caasf_is_const) fp_w += scatterer.caasf.c();
      fp_w *= w;
      if (have_fdp) fdp_w *= w;
      bool is_centric = space_group.is_centric();
      bool is_origin_centric = (
        is_centric ? space_group.is_origin_centric() : false);
      for(std::size_t i=0;i<miller_indices.size();i++) {
        miller::index<> const& h = miller_indices[i];
        f_t d_star_sq = unit_cell.d_star_sq(h);
        direct_one_h_one_scatterer<CosSinType, ScattererType> sf(
          cos_sin,
          space_group,
          h,
          d_star_sq,
          scatterer);
        if (is_origin_centric) {
          sf.f_calc = c_t(2*sf.f_calc.real(),0);
        }
        else if (is_centric) {
          f_t ht = f_t(h * space_group.inv_t()) / space_group.t_den();
          sf.f_calc += std::conj(sf.f_calc) * cos_sin.get(ht);
        }
        if (caasf_is_const) {
          if (have_fdp) sf.f_calc *= c_t(fp_w, fdp_w);
          else sf.f_calc *= fp_w;
        }
        else {
          f_t f0_w = scatterer.caasf.at_d_star_sq(d_star_sq) * w;
          if (have_fdp) sf.f_calc *= c_t(f0_w + fp_w, fdp_w);
          else sf.f_calc *= f0_w + fp_w;
        }
        f_calc[i] += sf.f_calc;
      }
    }
  };

  template <typename ScattererType=scatterer<> >
  class direct
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;

      direct() {}

      direct(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices, scatterers);
      }

      direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices, scatterers);
      }

      af::shared<std::complex<float_type> >
      f_calc() const { return f_calc_; }

    protected:
      af::shared<std::complex<float_type> > f_calc_;

      template <typename CosSinType>
      void
      compute(
        CosSinType const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers)
      {
        typedef float_type f_t;
        typedef std::complex<float_type> c_t;
        f_calc_.resize(miller_indices.size());
        af::ref<std::complex<float_type> > f_calc_ref = f_calc_.ref();
        for(std::size_t i=0;i<scatterers.size();i++) {
          ScattererType const& scatterer = scatterers[i];
          direct_one_scatterer<CosSinType, ScattererType> sf(
            cos_sin, unit_cell, space_group, miller_indices, scatterer,
            f_calc_ref);
        }
        if (space_group.n_ltr() > 1) {
          for(std::size_t i=0;i<miller_indices.size();i++) {
            f_calc_ref[i] *= space_group.n_ltr();
          }
        }
      }
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
