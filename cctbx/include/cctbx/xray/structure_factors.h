/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_XRAY_STRUCTURE_FACTORS_H
#define CCTBX_XRAY_STRUCTURE_FACTORS_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/sgtbx/miller_ops.h>

namespace cctbx { namespace xray { namespace structure_factors {

  /* d(f_calc)/d(site[i]) = two_pi * f_calc_hr * (j * hr[i])
     d(f_calc)/d(u_iso) = -two_pi_sq * f_calc_h * d_star_sq
     d(f_calc)/d(u_star(k,l)) = -two_pi_sq * f_calc_hr * (hr[k] * hr[l])
     d(f_calc)/d(occupancy) = f_calc_h / occupancy
     d(f_calc)/d(fp) = f_calc_h / (f0 + fp + j * fdp)
     d(f_calc)/d(fdp) = j f_calc_h / (f0 + fp + j * fdp)
   */

  template <typename ScattererType = scatterer<> >
  struct direct_with_first_derivatives_one_h_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    direct_with_first_derivatives_one_h_one_scatterer(
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      float_type d_star_sq,
      ScattererType const& scatterer,
      const std::complex<float_type>* d_target_d_f_calc,
      bool d_site_flag,
      bool d_u_star_flag)
    :
      const_h_sum(0,0)
    {
      typedef float_type f_t;
      typedef std::complex<f_t> c_t;
      f_t two_pi(scitbx::constants::two_pi);
      if (d_site_flag) d_target_d_site.fill(0);
      if (d_u_star_flag) d_target_d_u_star.fill(0);
      fractional<float_type> dtds_term;
      scitbx::sym_mat3<f_t> dw_coeff;
      f_t f0 = scatterer.caasf.at_d_star_sq(d_star_sq);
      f0_fp_fdp = f0 + scatterer.fp_fdp;
      f0_fp_fdp_w = f0_fp_fdp * scatterer.weight();
      for(std::size_t s=0;s<space_group.n_smx();s++) {
        miller::index<> hr = h * space_group.smx(s).r();
        f_t hrx = hr * scatterer.site;
        sgtbx::tr_vec t = space_group.smx(s).t();
        if (d_u_star_flag) {
          dw_coeff = adptbx::debye_waller_factor_u_star_coefficients(
            hr, scitbx::type_holder<f_t>());
        }
        c_t sum_inv(0,0);
        if (d_site_flag) dtds_term.fill(0);
        for(std::size_t i=0;i<space_group.f_inv();i++) {
          if (i) {
            hr = -hr;
            hrx = -hrx;
            t = space_group.inv_t() - t;
          }
          c_t sum_ltr(0,0);
          for(std::size_t l=0;l<space_group.n_ltr();l++) {
            f_t ht = f_t(h * (t + space_group.ltr(l))) / space_group.t_den();
            f_t phase = two_pi * (hrx + ht);
            c_t e_j_phase(std::cos(phase), std::sin(phase));
            sum_ltr += e_j_phase;
          }
          if (d_site_flag) {
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
          if (d_site_flag) dtds_term *= dw;
        }
        if (d_site_flag) d_target_d_site += dtds_term;
        if (d_u_star_flag) {
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
        if (d_site_flag) d_target_d_site *= dw;
      }
    }

    std::complex<float_type> f0_fp_fdp;
    std::complex<float_type> f0_fp_fdp_w;
    std::complex<float_type> const_h_sum;
    fractional<float_type> d_target_d_site;
    scitbx::sym_mat3<float_type> d_target_d_u_star;
  };

  template <typename ScattererType = scatterer<> >
  struct direct_with_first_derivatives_one_scatterer
  {
    typedef typename ScattererType::float_type float_type;

    direct_with_first_derivatives_one_scatterer(
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      af::const_ref<miller::index<> > const& miller_indices,
      ScattererType const& scatterer,
      af::ref<std::complex<float_type> > const& f_calc,
      af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
      bool d_site_flag,
      bool d_u_iso_flag,
      bool d_u_star_flag,
      bool d_occupancy_flag,
      bool d_fp_flag,
      bool d_fdp_flag)
    :
      d_target_d_site(0,0,0),
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
        direct_with_first_derivatives_one_h_one_scatterer<ScattererType> sf(
          space_group,
          h,
          d_star_sq,
          scatterer,
          d_t_d_f,
          d_site_flag,
          d_u_star_flag);
        f_calc[i] += sf.const_h_sum * sf.f0_fp_fdp_w;
        if (d_t_d_f) {
          if (d_site_flag) d_target_d_site += sf.d_target_d_site;
          if (d_u_star_flag) d_target_d_u_star += sf.d_target_d_u_star;
          if (d_u_iso_flag || d_occupancy_flag) {
            c_t t = sf.const_h_sum * sf.f0_fp_fdp;
            f_t d = d_t_d_f->real() * t.real()
                  + d_t_d_f->imag() * t.imag();
            d *= scatterer.weight_without_occupancy();
            if (d_u_iso_flag) {
              d_target_d_u_iso += d * scatterer.occupancy * d_star_sq;
            }
            if (d_occupancy_flag) {
              d_target_d_occupancy += d;
            }
          }
          if (d_fp_flag || d_fdp_flag) {
            c_t f = sf.const_h_sum * scatterer.weight();
            if (d_fp_flag) {
              d_target_d_fp += d_t_d_f->real() * f.real()
                             + d_t_d_f->imag() * f.imag();
            }
            if (d_fdp_flag) {
              d_target_d_fdp += d_t_d_f->imag() * f.real()
                              - d_t_d_f->real() * f.imag();
            }
          }
          d_t_d_f++;
        }
      }
      if (d_site_flag) d_target_d_site *= scitbx::constants::two_pi;
      if (d_u_iso_flag) d_target_d_u_iso *= -scitbx::constants::two_pi_sq;
      if (d_u_star_flag) d_target_d_u_star *= -scitbx::constants::two_pi_sq;
    }

    fractional<float_type> d_target_d_site;
    float_type d_target_d_u_iso;
    scitbx::sym_mat3<float_type> d_target_d_u_star;
    float_type d_target_d_occupancy;
    float_type d_target_d_fp;
    float_type d_target_d_fdp;
  };

  template <typename ScattererType = scatterer<> >
  class direct_with_first_derivatives
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;

      direct_with_first_derivatives() {}

      direct_with_first_derivatives(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        bool d_site_flag,
        bool d_u_iso_flag,
        bool d_u_star_flag,
        bool d_occupancy_flag,
        bool d_fp_flag,
        bool d_fdp_flag)
      {
        CCTBX_ASSERT(d_target_d_f_calc.size() == 0 ||
                     d_target_d_f_calc.size() == miller_indices.size());
        CCTBX_ASSERT(!(   d_site_flag
                       || d_u_iso_flag
                       || d_u_star_flag
                       || d_occupancy_flag
                       || d_fp_flag
                       || d_fdp_flag)
                  || d_target_d_f_calc.size() == miller_indices.size());
        f_calc_.resize(miller_indices.size());
        if (d_site_flag) d_target_d_site_.reserve(scatterers.size());
        if (d_u_iso_flag) d_target_d_u_iso_.reserve(scatterers.size());
        if (d_u_star_flag) d_target_d_u_star_.reserve(scatterers.size());
        if (d_occupancy_flag) d_target_d_occupancy_.reserve(scatterers.size());
        if (d_fp_flag) d_target_d_fp_.reserve(scatterers.size());
        if (d_fdp_flag) d_target_d_fdp_.reserve(scatterers.size());
        for(std::size_t i=0;i<scatterers.size();i++) {
          ScattererType const& scatterer = scatterers[i];
          direct_with_first_derivatives_one_scatterer<ScattererType> sf(
            unit_cell,
            space_group,
            miller_indices,
            scatterer,
            f_calc_.ref(),
            d_target_d_f_calc,
            d_site_flag,
            d_u_iso_flag && !scatterer.anisotropic_flag,
            d_u_star_flag && scatterer.anisotropic_flag,
            d_occupancy_flag,
            d_fp_flag,
            d_fdp_flag);
          if (d_site_flag) {
            d_target_d_site_.push_back(sf.d_target_d_site);
          }
          if (d_u_iso_flag) {
            d_target_d_u_iso_.push_back(sf.d_target_d_u_iso);
          }
          if (d_u_star_flag) {
            d_target_d_u_star_.push_back(sf.d_target_d_u_star);
          }
          if (d_occupancy_flag) {
            d_target_d_occupancy_.push_back(sf.d_target_d_occupancy);
          }
          if (d_fp_flag) {
            d_target_d_fp_.push_back(sf.d_target_d_fp);
          }
          if (d_fdp_flag) {
            d_target_d_fdp_.push_back(sf.d_target_d_fdp);
          }
        }
      }

      af::shared<std::complex<float_type> >
      f_calc() const { return f_calc_; }

      af::shared<scitbx::vec3<float_type> >
      d_target_d_site() const { return d_target_d_site_; }

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
      af::shared<std::complex<float_type> > f_calc_;
      af::shared<scitbx::vec3<float_type> > d_target_d_site_;
      af::shared<float_type> d_target_d_u_iso_;
      af::shared<scitbx::sym_mat3<float_type> > d_target_d_u_star_;
      af::shared<float_type> d_target_d_occupancy_;
      af::shared<float_type> d_target_d_fp_;
      af::shared<float_type> d_target_d_fdp_;
  };

  template <typename FloatType>
  void
  d_target_d_site_in_place_frac_as_cart(
    uctbx::unit_cell const& unit_cell,
    af::ref<scitbx::vec3<FloatType> > const& d_target_d_site)
  {
    for(std::size_t i=0;i<d_target_d_site.size();i++) {
      d_target_d_site[i] = d_target_d_site[i]
                         * unit_cell.fractionalization_matrix();
    }
  }

}}} // namespace cctbx::xray::structure_factors

#endif // CCTBX_XRAY_STRUCTURE_FACTORS_H
