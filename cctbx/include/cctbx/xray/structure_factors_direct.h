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
      for(std::size_t s=0;s<space_group.n_smx();s++) {
        miller::index<> hr = h * space_group.smx(s).r();
        f_t hrx = hr * scatterer.site;
        sgtbx::tr_vec t = space_group.smx(s).t();
        c_t sum_inv(0,0);
        for(std::size_t i=0;i<space_group.f_inv();i++) {
          if (i) {
            hr = -hr;
            hrx = -hrx;
            t = space_group.inv_t() - t;
          }
          for(std::size_t l=0;l<space_group.n_ltr();l++) {
            f_t ht = f_t(h * (t + space_group.ltr(l))) / space_group.t_den();
            sum_inv += cos_sin.get(hrx + ht);
          }
        }
        if (scatterer.anisotropic_flag) {
          f_t dw = adptbx::debye_waller_factor_u_star(hr, scatterer.u_star);
          sum_inv *= dw;
        }
        f_calc += sum_inv;
      }
      if (!scatterer.anisotropic_flag && scatterer.u_iso != 0) {
        f_t dw=adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        f_calc *= dw;
      }
      f_t f0 = scatterer.caasf.at_d_star_sq(d_star_sq);
      f_calc *= (f0 + scatterer.fp_fdp) * scatterer.weight();
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
      for(std::size_t i=0;i<miller_indices.size();i++) {
        miller::index<> const& h = miller_indices[i];
        f_t d_star_sq = unit_cell.d_star_sq(h);
        direct_one_h_one_scatterer<CosSinType, ScattererType> sf(
          cos_sin,
          space_group,
          h,
          d_star_sq,
          scatterer);
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
        f_calc_.resize(miller_indices.size());
        for(std::size_t i=0;i<scatterers.size();i++) {
          ScattererType const& scatterer = scatterers[i];
          direct_one_scatterer<math::cos_sin_exact<float_type>,
                               ScattererType> sf(
            cos_sin,
            unit_cell,
            space_group,
            miller_indices,
            scatterer,
            f_calc_.ref());
        }
      }

      direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers)
      {
        f_calc_.resize(miller_indices.size());
        for(std::size_t i=0;i<scatterers.size();i++) {
          ScattererType const& scatterer = scatterers[i];
          direct_one_scatterer<math::cos_sin_table<float_type>,
                               ScattererType> sf(
            cos_sin,
            unit_cell,
            space_group,
            miller_indices,
            scatterer,
            f_calc_.ref());
        }
      }

      af::shared<std::complex<float_type> >
      f_calc() const { return f_calc_; }

    protected:
      af::shared<std::complex<float_type> > f_calc_;
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
