#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H

#include <cctbx/xray/scattering_dictionary.h>
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
  struct direct_one_h
  {
    typedef typename ScattererType::float_type float_type;

    direct_one_h(
      CosSinType const& cos_sin,
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      af::const_ref<ScattererType> const& scatterers,
      scattering_dictionary const& scattering_dict)
    :
      f_calc(0,0)
    {
      typedef scattering_dictionary::dict_type dict_type;
      typedef dict_type::const_iterator dict_iter;
      typedef float_type f_t;
      typedef std::complex<float_type> c_t;
      bool is_centric = space_group.is_centric();
      bool is_origin_centric = (
        is_centric ? space_group.is_origin_centric() : false);
      f_t d_star_sq = unit_cell.d_star_sq(h);
      dict_type const& scd = scattering_dict.dict();
      for(dict_iter di=scd.begin();di!=scd.end();di++) {
        bool caasf_is_const = di->second.coefficients.n_ab() == 0;
        f_t caasf_c = di->second.coefficients.c();
        f_t f0 = di->second.coefficients.at_d_star_sq(d_star_sq);
        af::const_ref<std::size_t>
          member_indices = di->second.member_indices.const_ref();
        for(std::size_t mi=0;mi<member_indices.size();mi++) {
          ScattererType const& scatterer = scatterers[member_indices[mi]];
          f_t w = scatterer.weight();
          f_t fp_w = scatterer.fp;
          f_t fdp_w = scatterer.fdp;
          bool have_fdp = fdp_w != 0;
          if (caasf_is_const) fp_w += caasf_c;
          fp_w *= w;
          if (have_fdp) fdp_w *= w;
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
            f_t f0_w = f0 * w;
            if (have_fdp) sf.f_calc *= c_t(f0_w + fp_w, fdp_w);
            else sf.f_calc *= f0_w + fp_w;
          }
          f_calc += sf.f_calc;
        }
      }
    }

    std::complex<float_type> f_calc;
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
        af::const_ref<ScattererType> const& scatterers,
        scattering_dictionary const& scattering_dict)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_dict);
      }

      direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        scattering_dictionary const& scattering_dict)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_dict);
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
        af::const_ref<ScattererType> const& scatterers,
        scattering_dictionary const& scattering_dict)
      {
        typedef float_type f_t;
        typedef std::complex<float_type> c_t;
        f_t n_ltr = static_cast<f_t>(space_group.n_ltr());
        f_calc_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          f_calc_.push_back(
            direct_one_h<CosSinType, ScattererType>(
              cos_sin, unit_cell, space_group, miller_indices[i],
              scatterers, scattering_dict).f_calc
            * n_ltr);
        }
      }
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
