#ifndef CCTBX_XRAY_GRADIENTS_DIRECT_H
#define CCTBX_XRAY_GRADIENTS_DIRECT_H

#include <cctbx/xray/scattering_type_registry.h>
#include <cctbx/xray/hr_ht_cache.h>
#include <cctbx/xray/packing_order.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <cctbx/math/cos_sin_table.h>

namespace cctbx { namespace xray { namespace structure_factors {

  /* d(f_calc)/d(site[i]) = two_pi * f_calc_hr * (j * hr[i])
     d(f_calc)/d(u_iso) = -two_pi_sq * f_calc_h * d_star_sq
     d(f_calc)/d(u_star(k,l)) = -two_pi_sq * f_calc_hr * (hr[k] * hr[l])
     d(f_calc)/d(occupancy) = f_calc_h / occupancy
     d(f_calc)/d(fp) = f_calc_h / (f0 + fp + j * fdp)
     d(f_calc)/d(fdp) = j f_calc_h / (f0 + fp + j * fdp)
   */

  template <typename CosSinType, typename ScattererType>
  struct gradients_direct_one_h_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    std::complex<float_type> f0_fp_fdp;
    std::complex<float_type> const_h_sum;
    fractional<float_type> d_target_d_site_frac;
    scitbx::sym_mat3<float_type> d_target_d_u_star;

    gradients_direct_one_h_one_scatterer(
      CosSinType const& cos_sin,
      sgtbx::space_group const& space_group,
      hr_ht_cache<float_type> const& hr_ht,
      float_type d_star_sq,
      float_type f0,
      ScattererType const& scatterer,
      std::complex<float_type> const& d_target_d_f_calc)
    :
      f0_fp_fdp(f0 + std::complex<float_type>(scatterer.fp, scatterer.fdp)),
      const_h_sum(0,0),
      d_target_d_site_frac(0,0,0),
      d_target_d_u_star(0,0,0,0,0,0)
    {
      typedef float_type f_t;
      typedef std::complex<f_t> c_t;
      scitbx::sym_mat3<f_t> dw_coeff;
      c_t f0_fp_fdp_w = f0_fp_fdp * scatterer.weight();
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
          if (scatterer.flags.grad_site()) {
            c_t f = f0_fp_fdp_w * term;
            f_t c = d_target_d_f_calc.imag() * f.real()
                  - d_target_d_f_calc.real() * f.imag();
            if (i) c *= -1;
            for(std::size_t i=0;i<3;i++) {
              dtds_term[i] += g.hr[i] * c;
            }
          }
          sum_inv += term;
        }
        f_t dw = (scatterer.flags.use_u_iso() ? dw_iso : 1);
        if (scatterer.flags.use_u_aniso()) {
          dw *= adptbx::debye_waller_factor_u_star(g.hr, scatterer.u_star);
        }
        sum_inv *= dw;
        if (scatterer.flags.grad_site()) {
          dtds_term *= dw;
          d_target_d_site_frac += dtds_term;
        }
        if (scatterer.flags.grad_u_aniso() && scatterer.flags.use_u_aniso()) {
          c_t f = f0_fp_fdp_w * sum_inv;
          f_t c = d_target_d_f_calc.real() * f.real()
                + d_target_d_f_calc.imag() * f.imag();
          d_target_d_u_star += dw_coeff * c;
        }
        const_h_sum += sum_inv;
      }
    }
  };

  namespace detail {

    template <typename FloatType>
    struct gradient_refs
    {
      gradient_refs(
        af::ref<scitbx::vec3<FloatType> > site_,
        af::ref<FloatType> u_iso_,
        af::ref<scitbx::sym_mat3<FloatType> > u_star_,
        af::ref<FloatType> occupancy_,
        af::ref<FloatType> fp_,
        af::ref<FloatType> fdp_)
      :
        site(site_),
        u_iso(u_iso_),
        u_star(u_star_),
        occupancy(occupancy_),
        fp(fp_),
        fdp(fdp_)
      {}

      af::ref<scitbx::vec3<FloatType> > site;
      af::ref<FloatType> u_iso;
      af::ref<scitbx::sym_mat3<FloatType> > u_star;
      af::ref<FloatType> occupancy;
      af::ref<FloatType> fp;
      af::ref<FloatType> fdp;
    };

    template <typename ElementType, typename FactorType>
    inline
    void
    in_place_multiply(af::ref<ElementType> const& array, FactorType factor)
    {
      ElementType* p = array.begin();
      ElementType* e = array.end();
      while (p != e) (*p++) *= factor;
    }

  } // namespace detail

  template <typename CosSinType, typename ScattererType>
  struct gradients_direct_one_h
  {
    typedef typename ScattererType::float_type float_type;

    gradients_direct_one_h(
      CosSinType const& cos_sin,
      sgtbx::space_group const& space_group,
      af::const_ref<ScattererType> const& scatterers,
      af::const_ref<std::size_t> scattering_type_indices,
      miller::index<> h,
      float_type d_star_sq,
      af::const_ref<double> const& form_factors,
      std::complex<float_type> const& d_target_d_f_calc,
      detail::gradient_refs<float_type> gr_refs)
    {
      typedef float_type f_t;
      typedef std::complex<float_type> c_t;
      hr_ht_cache<f_t> hr_ht(space_group, h);
      for(std::size_t i_sc=0;i_sc<scatterers.size();i_sc++) {
        ScattererType const& scatterer = scatterers[i_sc];
        f_t f0 = form_factors[scattering_type_indices[i_sc]];
        bool gf_u_iso = scatterer.flags.grad_u_iso() &&
                        scatterer.flags.use_u_iso();
        bool gf_u_aniso = scatterer.flags.grad_u_aniso() &&
                          scatterer.flags.use_u_aniso();
        gradients_direct_one_h_one_scatterer<CosSinType, ScattererType> sf(
          cos_sin,
          space_group,
          hr_ht,
          d_star_sq,
          f0,
          scatterer,
          d_target_d_f_calc);
        if (scatterer.flags.grad_site()) {
          gr_refs.site[i_sc] += sf.d_target_d_site_frac;
        }
        if (gf_u_aniso) {
          gr_refs.u_star[i_sc] += sf.d_target_d_u_star;
        }
        if (gf_u_iso || scatterer.flags.grad_occupancy()) {
          c_t t = sf.const_h_sum * sf.f0_fp_fdp;
          f_t d = d_target_d_f_calc.real() * t.real()
                + d_target_d_f_calc.imag() * t.imag();
          d *= scatterer.weight_without_occupancy();
          if (gf_u_iso) {
            gr_refs.u_iso[i_sc] += d * scatterer.occupancy * d_star_sq;
          }
          if (scatterer.flags.grad_occupancy()) {
            gr_refs.occupancy[i_sc] += d;
          }
        }
        if (scatterer.flags.grad_fp() || scatterer.flags.grad_fdp()) {
          c_t f = sf.const_h_sum * scatterer.weight();
          if (scatterer.flags.grad_fp()) {
            gr_refs.fp[i_sc] += d_target_d_f_calc.real() * f.real()
                              + d_target_d_f_calc.imag() * f.imag();
          }
          if (scatterer.flags.grad_fdp()) {
            gr_refs.fdp[i_sc] += d_target_d_f_calc.imag() * f.real()
                               - d_target_d_f_calc.real() * f.imag();
          }
        }
      }
    }
  };

  template <typename ScattererType=scatterer<> >
  class gradients_direct
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;

      gradients_direct() {}

      gradients_direct(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<float_type> const& u_iso_refinable_params,
        xray::scattering_type_registry const& scattering_type_registry,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        std::size_t n_parameters=0)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, u_iso_refinable_params,
                scattering_type_registry, site_symmetry_table,
                d_target_d_f_calc, n_parameters);
      }

      gradients_direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<float_type> const& u_iso_refinable_params,
        xray::scattering_type_registry const& scattering_type_registry,
        sgtbx::site_symmetry_table const& site_symmetry_table,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        std::size_t n_parameters=0)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, u_iso_refinable_params,
                scattering_type_registry, site_symmetry_table,
                d_target_d_f_calc, n_parameters);
      }

      af::shared<float_type>
      packed() const { return packed_; }

      af::shared<scitbx::vec3<float_type> >
      d_target_d_site_frac() const { return d_target_d_site_frac_; }

      af::shared<float_type>
      d_target_d_u_iso() const { return d_target_d_u_iso_; }

      af::shared<scitbx::sym_mat3<float_type> >
      d_target_d_u_star() const { return d_target_d_u_star_; }

      af::shared<float_type>
      d_target_d_occupancy() const { return d_target_d_occupancy_; }

      af::shared<float_type>
      d_target_d_fp() const { return d_target_d_fp_; }

      af::shared<float_type>
      d_target_d_fdp() const { return d_target_d_fdp_; }

    protected:
      af::shared<float_type> packed_;
      af::shared<scitbx::vec3<float_type> > d_target_d_site_frac_;
      af::shared<float_type> d_target_d_u_iso_;
      af::shared<scitbx::sym_mat3<float_type> > d_target_d_u_star_;
      af::shared<float_type> d_target_d_occupancy_;
      af::shared<float_type> d_target_d_fp_;
      af::shared<float_type> d_target_d_fdp_;

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
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        std::size_t n_parameters)
      {
        CCTBX_ASSERT(d_target_d_f_calc.size() == miller_indices.size());
        cctbx::xray::scatterer_grad_flags_counts grad_flags_counts(scatterers);
        if(grad_flags_counts.tan_u_iso != 0 && grad_flags_counts.u_iso != 0) {
          CCTBX_ASSERT(u_iso_refinable_params.size() == scatterers.size());
        }
        CCTBX_ASSERT(grad_flags_counts.n_parameters() != 0);
        typedef float_type f_t;
        if (n_parameters != 0) {
          packed_.reserve(n_parameters);
        }
        if (grad_flags_counts.site != 0) {
          d_target_d_site_frac_.resize(
            scatterers.size(), scitbx::vec3<f_t>(0,0,0));
        }
        if (grad_flags_counts.u_iso != 0) {
          d_target_d_u_iso_.resize(scatterers.size(), 0);
        }
        if (grad_flags_counts.u_aniso != 0) {
          d_target_d_u_star_.resize(
            scatterers.size(), scitbx::sym_mat3<f_t>(0,0,0,0,0,0));
        }
        if (grad_flags_counts.occupancy != 0) {
          d_target_d_occupancy_.resize(scatterers.size(), 0);
        }
        if (grad_flags_counts.fp != 0) {
          d_target_d_fp_.resize(scatterers.size(), 0);
        }
        if (grad_flags_counts.fdp != 0) {
          d_target_d_fdp_.resize(scatterers.size(), 0);
        }
        detail::gradient_refs<f_t> gr_refs(
          d_target_d_site_frac_.ref(),
          d_target_d_u_iso_.ref(),
          d_target_d_u_star_.ref(),
          d_target_d_occupancy_.ref(),
          d_target_d_fp_.ref(),
          d_target_d_fdp_.ref());
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
          gradients_direct_one_h<CosSinType, ScattererType>(
            cos_sin, space_group, scatterers, scattering_type_indices,
            h, d_star_sq, form_factors,
            d_target_d_f_calc[i], gr_refs);
        }
        f_t n_ltr = static_cast<f_t>(space_group.n_ltr());
        if (grad_flags_counts.site != 0) {
          detail::in_place_multiply(gr_refs.site,
            n_ltr * static_cast<f_t>(scitbx::constants::two_pi));
        }
        if (grad_flags_counts.u_iso != 0) {
          detail::in_place_multiply(gr_refs.u_iso,
            n_ltr * static_cast<f_t>(-scitbx::constants::two_pi_sq));
        }
        if (grad_flags_counts.u_aniso != 0) {
          detail::in_place_multiply(gr_refs.u_star,
            n_ltr * static_cast<f_t>(-scitbx::constants::two_pi_sq));
        }
        if (grad_flags_counts.occupancy != 0) {
          detail::in_place_multiply(gr_refs.occupancy, n_ltr);
        }
        if (grad_flags_counts.fp != 0) {
          detail::in_place_multiply(gr_refs.fp, n_ltr);
        }
        if (grad_flags_counts.fdp != 0) {
          detail::in_place_multiply(gr_refs.fdp, n_ltr);
        }
        if (grad_flags_counts.u_iso != 0 && grad_flags_counts.tan_u_iso != 0) {
          float_type* d_t_d_u = d_target_d_u_iso_.begin();
          for(std::size_t i=0;i<scatterers.size();i++,d_t_d_u++) {
            ScattererType const& scatterer = scatterers[i];
            if (   scatterer.flags.tan_u_iso()
                && scatterer.flags.use_u_iso()
                && scatterer.flags.param > 0
                && scatterer.flags.grad_u_iso()) {
              f_t pi = scitbx::constants::pi;
              f_t u_iso_max = adptbx::b_as_u(
                                      static_cast<f_t>(scatterer.flags.param));
              (*d_t_d_u) *= u_iso_max/pi/(1.+u_iso_refinable_params[i]*
                                             u_iso_refinable_params[i]);
            }
          }
        }
        if (grad_flags_counts.site != 0) {
          average_special_position_site_gradients(
            site_symmetry_table,
            d_target_d_site_frac_.ref());
        }
        if (n_parameters != 0) {
          BOOST_STATIC_ASSERT(packing_order_convention == 2);
          for(std::size_t i=0;i<scatterers.size();i++) {
            ScattererType const& scatterer = scatterers[i];
            if (scatterer.flags.grad_site()) {
              scitbx::vec3<f_t> d_target_d_site_cart =
                unit_cell.v_times_fractionalization_matrix_transpose(
                  /* v */ gr_refs.site[i]);
              for(std::size_t j=0;j<3;j++) {
                packed_.push_back(d_target_d_site_cart[j]);
              }
            }
            if(scatterer.flags.use_u_iso()) {
              if (scatterer.flags.grad_u_iso()) {
                packed_.push_back(gr_refs.u_iso[i]);
              }
            }
            if(scatterer.flags.use_u_aniso()) {
               if(scatterer.flags.grad_u_aniso()) {
                 scitbx::sym_mat3<double> d_target_d_u_cart =
                   adptbx::grad_u_star_as_u_cart(
                     unit_cell, gr_refs.u_star[i]);
                 for(std::size_t j=0;j<6;j++) {
                   packed_.push_back(d_target_d_u_cart[j]);
                 }
               }
             }
            if (scatterer.flags.grad_occupancy()) {
              packed_.push_back(gr_refs.occupancy[i]);
            }
            // XXX use scatterer flags and not grad_flags_counts.
            // XXX what's below may not work if selected part is refining
            if (grad_flags_counts.fp != 0) {
              packed_.push_back(gr_refs.fp[i]);
            }
            if (grad_flags_counts.fdp != 0) {
              packed_.push_back(gr_refs.fdp[i]);
            }
          }
          CCTBX_ASSERT(packed_.size() == n_parameters);
          d_target_d_site_frac_ = af::shared<scitbx::vec3<f_t> >();
          d_target_d_u_iso_ = af::shared<f_t>();
          d_target_d_u_star_ = af::shared<scitbx::sym_mat3<f_t> >();
          d_target_d_occupancy_ = af::shared<f_t>();
          d_target_d_fp_ = af::shared<f_t>();
          d_target_d_fdp_ = af::shared<f_t>();
        }
      }
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_GRADIENTS_DIRECT_H
