#ifndef CCTBX_XRAY_GRADIENTS_DIRECT_H
#define CCTBX_XRAY_GRADIENTS_DIRECT_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/xray/gradient_flags.h>
#include <cctbx/xray/packing_order.h>
#include <cctbx/math/cos_sin_table.h>
#include <cctbx/sgtbx/miller_ops.h>

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

    gradients_direct_one_h_one_scatterer(
      CosSinType const& cos_sin,
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      float_type d_star_sq,
      ScattererType const& scatterer,
      const std::complex<float_type>* d_target_d_f_calc,
      bool grad_flags_site,
      bool grad_flags_u_aniso)
    :
      const_h_sum(0,0)
    {
      typedef float_type f_t;
      typedef std::complex<f_t> c_t;
      if (grad_flags_site) d_target_d_site_frac.fill(0);
      if (grad_flags_u_aniso) d_target_d_u_star.fill(0);
      fractional<float_type> dtds_term;
      scitbx::sym_mat3<f_t> dw_coeff;
      f_t f0 = scatterer.caasf.at_d_star_sq(d_star_sq);
      f0_fp_fdp = f0 + scatterer.fp_fdp;
      f0_fp_fdp_w = f0_fp_fdp * scatterer.weight();
      for(std::size_t s=0;s<space_group.n_smx();s++) {
        miller::index<> hr = h * space_group.smx(s).r();
        f_t hrx = hr * scatterer.site;
        sgtbx::tr_vec t = space_group.smx(s).t();
        if (grad_flags_u_aniso) {
          dw_coeff = adptbx::debye_waller_factor_u_star_coefficients(
            hr, scitbx::type_holder<f_t>());
        }
        c_t sum_inv(0,0);
        if (grad_flags_site) dtds_term.fill(0);
        for(std::size_t i=0;i<space_group.f_inv();i++) {
          if (i) {
            hr = -hr;
            hrx = -hrx;
            t = space_group.inv_t() - t;
          }
          c_t sum_ltr(0,0);
          for(std::size_t l=0;l<space_group.n_ltr();l++) {
            f_t ht = f_t(h * (t + space_group.ltr(l))) / space_group.t_den();
            sum_ltr += cos_sin.get(hrx + ht);
          }
          if (grad_flags_site) {
            c_t f = f0_fp_fdp_w * sum_ltr;
            f_t c = d_target_d_f_calc->imag() * f.real()
                  - d_target_d_f_calc->real() * f.imag();
            for(std::size_t i=0;i<3;i++) {
              dtds_term[i] += hr[i] * c;
            }
          }
          sum_inv += sum_ltr;
        }
        if (scatterer.anisotropic_flag) {
          f_t dw = adptbx::debye_waller_factor_u_star(hr, scatterer.u_star);
          sum_inv *= dw;
          if (grad_flags_site) dtds_term *= dw;
        }
        if (grad_flags_site) d_target_d_site_frac += dtds_term;
        if (grad_flags_u_aniso) {
          c_t f = f0_fp_fdp_w * sum_inv;
          f_t c = d_target_d_f_calc->real() * f.real()
                + d_target_d_f_calc->imag() * f.imag();
          d_target_d_u_star += dw_coeff * c;
        }
        const_h_sum += sum_inv;
      }
      if (!scatterer.anisotropic_flag && scatterer.u_iso != 0) {
        f_t dw=adptbx::debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso);
        const_h_sum *= dw;
        if (grad_flags_site) d_target_d_site_frac *= dw;
      }
    }

    std::complex<float_type> f0_fp_fdp;
    std::complex<float_type> f0_fp_fdp_w;
    std::complex<float_type> const_h_sum;
    fractional<float_type> d_target_d_site_frac;
    scitbx::sym_mat3<float_type> d_target_d_u_star;
  };

  template <typename CosSinType, typename ScattererType>
  struct gradients_direct_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    gradients_direct_one_scatterer(
      CosSinType const& cos_sin,
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      af::const_ref<miller::index<> > const& miller_indices,
      ScattererType const& scatterer,
      af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
      gradient_flags const& grad_flags)
    :
      d_target_d_site_frac(0,0,0),
      d_target_d_u_iso(0),
      d_target_d_u_star(0,0,0,0,0,0),
      d_target_d_occupancy(0),
      d_target_d_fp(0),
      d_target_d_fdp(0)
    {
      typedef float_type f_t;
      typedef std::complex<float_type> c_t;
      const c_t* d_t_d_f = d_target_d_f_calc.begin();
      for(std::size_t i=0;i<miller_indices.size();i++) {
        miller::index<> const& h = miller_indices[i];
        f_t d_star_sq = unit_cell.d_star_sq(h);
        gradients_direct_one_h_one_scatterer<CosSinType, ScattererType> sf(
          cos_sin,
          space_group,
          h,
          d_star_sq,
          scatterer,
          d_t_d_f,
          grad_flags.site,
          grad_flags.u_aniso);
        if (d_t_d_f) {
          if (grad_flags.site) d_target_d_site_frac += sf.d_target_d_site_frac;
          if (grad_flags.u_aniso) d_target_d_u_star += sf.d_target_d_u_star;
          if (grad_flags.u_iso || grad_flags.occupancy) {
            c_t t = sf.const_h_sum * sf.f0_fp_fdp;
            f_t d = d_t_d_f->real() * t.real()
                  + d_t_d_f->imag() * t.imag();
            d *= scatterer.weight_without_occupancy();
            if (grad_flags.u_iso) {
              d_target_d_u_iso += d * scatterer.occupancy * d_star_sq;
            }
            if (grad_flags.occupancy) {
              d_target_d_occupancy += d;
            }
          }
          if (grad_flags.fp || grad_flags.fdp) {
            c_t f = sf.const_h_sum * scatterer.weight();
            if (grad_flags.fp) {
              d_target_d_fp += d_t_d_f->real() * f.real()
                             + d_t_d_f->imag() * f.imag();
            }
            if (grad_flags.fdp) {
              d_target_d_fdp += d_t_d_f->imag() * f.real()
                              - d_t_d_f->real() * f.imag();
            }
          }
          d_t_d_f++;
        }
      }
      if (grad_flags.site) d_target_d_site_frac *= scitbx::constants::two_pi;
      if (grad_flags.u_iso) d_target_d_u_iso *= -scitbx::constants::two_pi_sq;
      if (grad_flags.u_aniso) d_target_d_u_star*=-scitbx::constants::two_pi_sq;
    }

    fractional<float_type> d_target_d_site_frac;
    float_type d_target_d_u_iso;
    scitbx::sym_mat3<float_type> d_target_d_u_star;
    float_type d_target_d_occupancy;
    float_type d_target_d_fp;
    float_type d_target_d_fdp;
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
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        gradient_flags const& grad_flags,
        std::size_t n_parameters=0)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices, scatterers,
                d_target_d_f_calc, grad_flags, n_parameters);
      }

      gradients_direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        gradient_flags const& grad_flags,
        std::size_t n_parameters=0)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices, scatterers,
                d_target_d_f_calc, grad_flags, n_parameters);
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

      template <typename CosSinType>
      void
      compute(
        CosSinType const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        gradient_flags const& grad_flags,
        std::size_t n_parameters)
      {
        CCTBX_ASSERT(   d_target_d_f_calc.size() == 0
                     || d_target_d_f_calc.size() == miller_indices.size());
        CCTBX_ASSERT(   grad_flags.all_false()
                     || d_target_d_f_calc.size() == miller_indices.size());
        if (n_parameters != 0) {
          packed_.reserve(n_parameters);
        }
        else {
          if (grad_flags.site)d_target_d_site_frac_.reserve(scatterers.size());
          if (grad_flags.u_iso) d_target_d_u_iso_.reserve(scatterers.size());
          if (grad_flags.u_aniso)d_target_d_u_star_.reserve(scatterers.size());
          if (grad_flags.occupancy) d_target_d_occupancy_.reserve(
                                      scatterers.size());
          if (grad_flags.fp) d_target_d_fp_.reserve(scatterers.size());
          if (grad_flags.fdp) d_target_d_fdp_.reserve(scatterers.size());
        }
        for(std::size_t i=0;i<scatterers.size();i++) {
          ScattererType const& scatterer = scatterers[i];
          gradients_direct_one_scatterer<CosSinType, ScattererType> sf(
            cos_sin,
            unit_cell,
            space_group,
            miller_indices,
            scatterer,
            d_target_d_f_calc,
            grad_flags.adjust(scatterer.anisotropic_flag));
          if (n_parameters != 0) {
            BOOST_STATIC_ASSERT(packing_order_convention == 1);
            if (grad_flags.site) {
              scitbx::vec3<float_type> d_target_d_site_cart =
                sf.d_target_d_site_frac * unit_cell.fractionalization_matrix();
              for(std::size_t i=0;i<3;i++) {
                packed_.push_back(d_target_d_site_cart[i]);
              }
            }
            if (!scatterer.anisotropic_flag) {
              if (grad_flags.u_iso) {
                packed_.push_back(sf.d_target_d_u_iso);
              }
            }
            else {
              if (grad_flags.u_aniso) {
                scitbx::sym_mat3<double> d_target_d_u_cart =
                  adptbx::grad_u_star_as_u_cart(
                    unit_cell, sf.d_target_d_u_star);
                for(std::size_t i=0;i<6;i++) {
                  packed_.push_back(d_target_d_u_cart[i]);
                }
              }
            }
            if (grad_flags.occupancy) {
              packed_.push_back(sf.d_target_d_occupancy);
            }
            if (grad_flags.fp) {
              packed_.push_back(sf.d_target_d_fp);
            }
            if (grad_flags.fdp) {
              packed_.push_back(sf.d_target_d_fdp);
            }
          }
          else {
            if (grad_flags.site) {
              d_target_d_site_frac_.push_back(sf.d_target_d_site_frac);
            }
            if (grad_flags.u_iso) {
              d_target_d_u_iso_.push_back(sf.d_target_d_u_iso);
            }
            if (grad_flags.u_aniso) {
              d_target_d_u_star_.push_back(sf.d_target_d_u_star);
            }
            if (grad_flags.occupancy) {
              d_target_d_occupancy_.push_back(sf.d_target_d_occupancy);
            }
            if (grad_flags.fp) {
              d_target_d_fp_.push_back(sf.d_target_d_fp);
            }
            if (grad_flags.fdp) {
              d_target_d_fdp_.push_back(sf.d_target_d_fdp);
            }
          }
        }
        if (n_parameters != 0) {
          CCTBX_ASSERT(packed_.size() == n_parameters);
        }
      }
  };

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_GRADIENTS_DIRECT_H
