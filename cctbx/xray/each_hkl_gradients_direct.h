#ifndef CCTBX_XRAY_EACH_HKL_GRADIENTS_DIRECT_H
#define CCTBX_XRAY_EACH_HKL_GRADIENTS_DIRECT_H

#include <cctbx/xray/gradients_direct.h>

namespace cctbx { namespace xray { namespace structure_factors {

  /* d(f_calc)/d(site[i]) = two_pi * f_calc_hr * (j * hr[i])
     d(f_calc)/d(u_iso) = -two_pi_sq * f_calc_h * d_star_sq
     d(f_calc)/d(u_star(k,l)) = -two_pi_sq * f_calc_hr * (hr[k] * hr[l])
     d(f_calc)/d(occupancy) = f_calc_h / occupancy
     d(f_calc)/d(fp) = f_calc_h / (f0 + fp + j * fdp)
     d(f_calc)/d(fdp) = j f_calc_h / (f0 + fp + j * fdp)
   */

  template <typename CosSinType, typename ScattererType>
  struct each_hkl_gradients_direct_one_h_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    std::complex<float_type> f0_fp_fdp;
    std::complex<float_type> const_h_sum;
    fractional<float_type> d_target_d_site_frac;
    scitbx::sym_mat3<float_type> d_target_d_u_star;

    each_hkl_gradients_direct_one_h_one_scatterer(
      CosSinType const& cos_sin,
      sgtbx::space_group const& space_group,
      hr_ht_cache<float_type> const& hr_ht,
      float_type d_star_sq,
      float_type f0,
      ScattererType const& scatterer)
    :
      f0_fp_fdp(f0 + std::complex<float_type>(scatterer.fp, scatterer.fdp)),
      const_h_sum(0,0),
      d_target_d_site_frac(0,0,0),
      d_target_d_u_star(0,0,0,0,0,0)
    {
      typedef float_type f_t;
      typedef std::complex<f_t> c_t;
      scitbx::sym_mat3<f_t> dw_coeff;
      f_t dw_iso = 0;
      if (scatterer.flags.use_u_iso()) {
        dw_iso = adptbx::debye_waller_factor_u_iso(
          d_star_sq/4, scatterer.u_iso);
      }
      for(std::size_t i=0;i<hr_ht.groups.size();i++) {
        hr_ht_group<f_t> const& g = hr_ht.groups[i];
        f_t hrx = g.hr * scatterer.site;
        f_t ht = g.ht;
        if (scatterer.flags.use_u_aniso()) {
          dw_coeff = adptbx::debye_waller_factor_u_star_gradient_coefficients(
            g.hr, scitbx::type_holder<f_t>());
        }
        c_t sum_inv(0,0);
        fractional<float_type> dtds_term(0,0,0);
        for(std::size_t i=0;i<space_group.f_inv();i++) {
          if (i) {
            hrx *= -1;
            ht = hr_ht.h_inv_t - ht;
          }
          c_t term = cos_sin.get(hrx + ht);
          sum_inv += term;
        }
        f_t dw = (scatterer.flags.use_u_iso() ? dw_iso : 1);
        if (scatterer.flags.use_u_aniso()) {
          dw *= adptbx::debye_waller_factor_u_star(g.hr, scatterer.u_star);
        }
        sum_inv *= dw;
        const_h_sum += sum_inv;
      }
    }
  };

  template <typename CosSinType, typename ScattererType>
  struct each_hkl_gradients_direct_one_h
  {
    typedef typename ScattererType::float_type float_type;
    typedef std::complex<float_type> c_t;
    c_t fp, fdp;

    each_hkl_gradients_direct_one_h(
      CosSinType const& cos_sin,
      sgtbx::space_group const& space_group,
      af::const_ref<ScattererType> const& scatterers,
      af::const_ref<std::size_t> scattering_type_indices,
      miller::index<> h,
      float_type d_star_sq,
      af::const_ref<double> const& form_factors)
    {
      typedef float_type f_t;
      typedef std::complex<float_type> c_t;
      hr_ht_cache<f_t> hr_ht(space_group, h);
      for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
        ScattererType const& scatterer = scatterers[i_sc];
        f_t f0 = form_factors[scattering_type_indices[i_sc]];
        each_hkl_gradients_direct_one_h_one_scatterer<CosSinType, ScattererType> sf(
          cos_sin,
          space_group,
          hr_ht,
          d_star_sq,
          f0,
          scatterer);
        if (scatterer.flags.grad_fp() || scatterer.flags.grad_fdp()) {
          c_t f = sf.const_h_sum * scatterer.weight();
          if (scatterer.flags.grad_fp()) {
            fp += f;
          }
          if (scatterer.flags.grad_fdp()) {
            fdp += c_t(0,1) * f;
          }
        }
      }
    }
  };

  template <typename ScattererType=scatterer<> >
  class each_hkl_gradients_direct
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;
      typedef std::complex<float_type> c_t;

      each_hkl_gradients_direct() {}

      each_hkl_gradients_direct(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<float_type> const& u_iso_refinable_params,
        xray::scattering_type_registry const& scattering_type_registry,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        std::size_t n_parameters=0)
      {
        SCITBX_ASSERT(scattering_type_registry.size()<=1);
        // The purpose of this class is to compute the Miller arrays giving
        // d(F_calc[H])/d(parameter), where parameter is e.g. fprime or fdoubleprime.
        // Therefore, it only makes sense to consider scatterer lists where all the
        // sites are considered to be of the same chemical type, such as all methionine sulfurs, etc.
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, u_iso_refinable_params,
                scattering_type_registry, site_symmetry_table,
                n_parameters);
      }

      each_hkl_gradients_direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<float_type> const& u_iso_refinable_params,
        xray::scattering_type_registry const& scattering_type_registry,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        std::size_t n_parameters=0)
      {
        SCITBX_ASSERT(scattering_type_registry.size()<=1);
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, u_iso_refinable_params,
                scattering_type_registry, site_symmetry_table,
                n_parameters);
      }

      af::shared<c_t >
      d_fcalc_d_fp() const { return d_fcalc_d_fp_; }

      af::shared<c_t >
      d_fcalc_d_fdp() const { return d_fcalc_d_fdp_; }

    protected:
      af::shared<std::complex<float_type> > d_fcalc_d_fp_;
      af::shared<std::complex<float_type> > d_fcalc_d_fdp_;

      // compensates for rounding errors
      static
      void
      average_special_position_site_gradients(
        sgtbx::site_symmetry_table const& site_symmetry_table,
        af::ref<scitbx::vec3<float_type> > gradients)
      {
        CCTBX_ASSERT(gradients.size()
                  == site_symmetry_table.indices_const_ref().size());
        af::const_ref<std::size_t> sp_indices = site_symmetry_table
          .special_position_indices().const_ref();
        for(std::size_t i_sp=0;i_sp<sp_indices.size();i_sp++) {
          std::size_t i_seq = sp_indices[i_sp];
          gradients[i_seq] = gradients[i_seq]
                           * site_symmetry_table.get(i_seq).special_op().r();
        }
      }

      template <typename CosSinType>
      void
      compute(
        CosSinType const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<float_type> const& u_iso_refinable_params,
        xray::scattering_type_registry const& scattering_type_registry,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        std::size_t n_parameters)
      {
        cctbx::xray::scatterer_grad_flags_counts grad_flags_counts(scatterers);
        if(grad_flags_counts.tan_u_iso != 0 && grad_flags_counts.u_iso != 0) {
          CCTBX_ASSERT(u_iso_refinable_params.size() == scatterers.size());
        }
        CCTBX_ASSERT(grad_flags_counts.n_parameters() != 0);
        typedef float_type f_t;
        if (grad_flags_counts.fp != 0) {
          d_fcalc_d_fp_.resize(miller_indices.size(), 0);
        }
        if (grad_flags_counts.fdp != 0) {
          d_fcalc_d_fdp_.resize(miller_indices.size(), 0);
        }
        af::shared<std::size_t>
          scattering_type_indices_memory
            = scattering_type_registry.unique_indices(scatterers);
        af::const_ref<std::size_t>
          scattering_type_indices = scattering_type_indices_memory.const_ref();
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller::index<> h = miller_indices[i];
          f_t d_star_sq = unit_cell.d_star_sq(h);
          af::shared<double>
            form_factors_memory
              = scattering_type_registry.unique_form_factors_at_d_star_sq(
                  d_star_sq);
          af::const_ref<double> form_factors = form_factors_memory.const_ref();
          each_hkl_gradients_direct_one_h<CosSinType, ScattererType> ONE_H(
            cos_sin, space_group, scatterers, scattering_type_indices,
            h, d_star_sq, form_factors);
          d_fcalc_d_fp_[i]=ONE_H.fp;
          d_fcalc_d_fdp_[i]=ONE_H.fdp;
        }
        f_t n_ltr = static_cast<f_t>(space_group.n_ltr());

        if (grad_flags_counts.fp != 0) {
          detail::in_place_multiply(d_fcalc_d_fp_.ref(), n_ltr);
        }
        if (grad_flags_counts.fdp != 0) {
          detail::in_place_multiply(d_fcalc_d_fdp_.ref(), n_ltr);
        }
      }
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_EACH_HKL_GRADIENTS_DIRECT_H
