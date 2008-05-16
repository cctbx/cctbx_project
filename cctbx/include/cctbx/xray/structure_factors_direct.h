#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H

#include <cctbx/xray/scattering_type_registry.h>
#include <cctbx/xray/hr_ht_cache.h>
#include <cctbx/math/cos_sin_table.h>

namespace cctbx { namespace xray { namespace structure_factors {

  template <typename CosSinType, typename ScattererType>
  struct direct_sum_over_equivalent_h
  {
    typedef typename ScattererType::float_type float_type;
    typedef std::complex<float_type> complex_type;

    direct_sum_over_equivalent_h(
      CosSinType const& cos_sin_,
      sgtbx::space_group const& space_group_,
      miller::index<> h,
      float_type d_star_sq_)
    :
      cos_sin(cos_sin_),
      hr_ht(cos_sin_, space_group_, h),
      d_star_sq(d_star_sq_),
      sum_f_calc(0,0)
    {}

    void add_contribution_of(ScattererType const& scatterer,
                             float_type f0)
    {
      typedef float_type f_t;
      typedef complex_type c_t;
      c_t f_calc(0,0);
      for(std::size_t i=0;i<hr_ht.groups.size();i++) {
        hr_ht_group<f_t> const& g = hr_ht.groups[i];
        f_t hrx = g.hr * scatterer.site;
        c_t term = cos_sin.get(hrx + g.ht);
        if (scatterer.flags.use_u_aniso()) {
          f_t dw = adptbx::debye_waller_factor_u_star(g.hr, scatterer.u_star);
          term *= dw;
        }
        f_calc += term;
      }
      if (hr_ht.is_origin_centric) {
        f_calc = c_t(2*f_calc.real(),0);
      }
      else if (hr_ht.is_centric) {
        f_calc += std::conj(f_calc) * hr_ht.f_h_inv_t;
      }
      if (scatterer.flags.use_u_iso() && scatterer.u_iso != 0) {
        f_t dw=adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        f_calc *= dw;
      }
      f_t w = scatterer.weight();
      f_t f0p_w = (f0 + scatterer.fp) * w;
      f_t fdp_w = scatterer.fdp;
      if (fdp_w != 0) {
        fdp_w *= w;
        f_calc *= c_t(f0p_w, fdp_w);
      }
      else {
        f_calc *= f0p_w;
      }
      sum_f_calc += f_calc;
    }

    complex_type f_calc() {
      return sum_f_calc * hr_ht.ltr_factor;
    }

    CosSinType const &cos_sin;
    hr_ht_cache<float_type> hr_ht;
    float_type d_star_sq;
    complex_type sum_f_calc;
  };


  template <class ScattererType=scatterer<>,
            class FormFactorType=eltbx::xray_scattering::gaussian>
  class direct
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;
      typedef xray::generic_scattering_type_registry<FormFactorType>
              scattering_type_registry_t;

      direct() {}

      direct(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        scattering_type_registry_t const& scattering_type_registry)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_type_registry);
      }

      direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        scattering_type_registry_t const& scattering_type_registry)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_type_registry);
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
        scattering_type_registry_t const& scattering_type_registry)
      {
        typedef float_type f_t;
        typedef std::complex<float_type> c_t;
        f_calc_.reserve(miller_indices.size());
        af::shared<std::size_t> scattering_type_indices
          = scattering_type_registry.unique_indices(scatterers);
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::index<> h = miller_indices[i];
          f_t d_star_sq = unit_cell.d_star_sq(h);
          af::shared<double> form_factors
            = scattering_type_registry.unique_form_factors_at_d_star_sq(
                d_star_sq);
          direct_sum_over_equivalent_h<CosSinType, ScattererType>
            sum(cos_sin, space_group, h, d_star_sq);
          for(std::size_t j=0; j<scatterers.size(); ++j) {
            sum.add_contribution_of(scatterers[j],
                                    form_factors[scattering_type_indices[j]]);
          }
          f_calc_.push_back(sum.f_calc());
        }
      }
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_STRUCTURE_FACTORS_DIRECT_H
