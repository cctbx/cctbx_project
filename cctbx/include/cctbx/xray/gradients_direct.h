#ifndef CCTBX_XRAY_GRADIENTS_DIRECT_H
#define CCTBX_XRAY_GRADIENTS_DIRECT_H

#include <cctbx/xray/scattering_dictionary.h>
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
      float_type f0,
      ScattererType const& scatterer,
      std::complex<float_type> const& d_target_d_f_calc,
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
      f0_fp_fdp = f0 + c_t(scatterer.fp, scatterer.fdp);
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
            f_t c = d_target_d_f_calc.imag() * f.real()
                  - d_target_d_f_calc.real() * f.imag();
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
          f_t c = d_target_d_f_calc.real() * f.real()
                + d_target_d_f_calc.imag() * f.imag();
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

  } // namespace detail

  template <typename CosSinType, typename ScattererType>
  struct gradients_direct_one_h
  {
    typedef typename ScattererType::float_type float_type;

    gradients_direct_one_h(
      CosSinType const& cos_sin,
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group,
      miller::index<> const& h,
      af::const_ref<ScattererType> const& scatterers,
      scattering_dictionary const& scattering_dict,
      std::complex<float_type> const& d_target_d_f_calc,
      gradient_flags const& grad_flags,
      detail::gradient_refs<float_type> gr_refs)
    {
      typedef scattering_dictionary::dict_type dict_type;
      typedef dict_type::const_iterator dict_iter;
      typedef float_type f_t;
      typedef std::complex<float_type> c_t;
      f_t d_star_sq = unit_cell.d_star_sq(h);
      dict_type const& scd = scattering_dict.dict();
      std::size_t gr_pos = 0;
      for(dict_iter di=scd.begin();di!=scd.end();di++) {
        f_t f0 = di->second.coefficients.at_d_star_sq(d_star_sq);
        af::const_ref<std::size_t>
          member_indices = di->second.member_indices.const_ref();
        for(std::size_t mi=0;mi<member_indices.size();mi++,gr_pos++) {
          ScattererType const& scatterer = scatterers[member_indices[mi]];
          bool gf_u_iso = grad_flags.u_iso && !scatterer.anisotropic_flag;
          bool gf_u_aniso = grad_flags.u_aniso && scatterer.anisotropic_flag;
          gradients_direct_one_h_one_scatterer<CosSinType, ScattererType> sf(
            cos_sin,
            space_group,
            h,
            d_star_sq,
            f0,
            scatterer,
            d_target_d_f_calc,
            grad_flags.site,
            gf_u_aniso);
          if (grad_flags.site) {
            gr_refs.site[gr_pos] += sf.d_target_d_site_frac;
          }
          if (gf_u_aniso) {
            gr_refs.u_star[gr_pos] += sf.d_target_d_u_star;
          }
          if (gf_u_iso || grad_flags.occupancy) {
            c_t t = sf.const_h_sum * sf.f0_fp_fdp;
            f_t d = d_target_d_f_calc.real() * t.real()
                  + d_target_d_f_calc.imag() * t.imag();
            d *= scatterer.weight_without_occupancy();
            if (gf_u_iso) {
              gr_refs.u_iso[gr_pos] += d * scatterer.occupancy * d_star_sq;
            }
            if (grad_flags.occupancy) {
              gr_refs.occupancy[gr_pos] += d;
            }
          }
          if (grad_flags.fp || grad_flags.fdp) {
            c_t f = sf.const_h_sum * scatterer.weight();
            if (grad_flags.fp) {
              gr_refs.fp[gr_pos] += d_target_d_f_calc.real() * f.real()
                                  + d_target_d_f_calc.imag() * f.imag();
            }
            if (grad_flags.fdp) {
              gr_refs.fdp[gr_pos] += d_target_d_f_calc.imag() * f.real()
                                   - d_target_d_f_calc.real() * f.imag();
            }
          }
        }
      }
    }
  };

  namespace detail {

    template <typename ElementType>
    af::shared<ElementType>
    unscramble(
      af::const_ref<ElementType> const& data,
      af::const_ref<std::size_t> const& permutation)
    {
      CCTBX_ASSERT(data.size() == permutation.size());
      af::shared<ElementType> result_mem;
      if (data.size() > 0) result_mem.resize(data.size(), data[0]);
      af::ref<ElementType> result = result_mem.ref();
      for(std::size_t i=0;i<data.size();i++) {
        CCTBX_ASSERT(permutation[i] < result.size());
        result[permutation[i]] = data[i];
      }
      return result_mem;
    }

    template <typename ElementType,
              typename FactorType>
    af::shared<ElementType>
    unscramble(
      af::const_ref<ElementType> const& data,
      af::const_ref<std::size_t> const& permutation,
      FactorType const& factor)
    {
      CCTBX_ASSERT(data.size() == permutation.size());
      af::shared<ElementType> result_mem;
      if (data.size() > 0) result_mem.resize(data.size(), data[0]);
      af::ref<ElementType> result = result_mem.ref();
      for(std::size_t i=0;i<data.size();i++) {
        CCTBX_ASSERT(permutation[i] < result.size());
        result[permutation[i]] = data[i] * factor;
      }
      return result_mem;
    }

  } // namespace detail

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
        scattering_dictionary const& scattering_dict,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        gradient_flags const& grad_flags,
        std::size_t n_parameters=0)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_dict,
                d_target_d_f_calc, grad_flags, n_parameters);
      }

      gradients_direct(
        math::cos_sin_table<float_type> const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        scattering_dictionary const& scattering_dict,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        gradient_flags const& grad_flags,
        std::size_t n_parameters=0)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_dict,
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
        scattering_dictionary const& scattering_dict,
        af::const_ref<std::complex<float_type> > const& d_target_d_f_calc,
        gradient_flags const& grad_flags,
        std::size_t n_parameters)
      {
        CCTBX_ASSERT(d_target_d_f_calc.size() == miller_indices.size());
        CCTBX_ASSERT(!grad_flags.all_false());
        typedef float_type f_t;
        if (n_parameters != 0) {
          packed_.reserve(n_parameters);
        }
        if (grad_flags.site) {
          d_target_d_site_frac_.resize(
            scatterers.size(), scitbx::vec3<f_t>(0,0,0));
        }
        if (grad_flags.u_iso) {
          d_target_d_u_iso_.resize(scatterers.size(), 0);
        }
        if (grad_flags.u_aniso) {
          d_target_d_u_star_.resize(
            scatterers.size(), scitbx::sym_mat3<f_t>(0,0,0,0,0,0));
        }
        if (grad_flags.occupancy) {
          d_target_d_occupancy_.resize(scatterers.size(), 0);
        }
        if (grad_flags.fp) {
          d_target_d_fp_.resize(scatterers.size(), 0);
        }
        if (grad_flags.fdp) {
          d_target_d_fdp_.resize(scatterers.size(), 0);
        }
        detail::gradient_refs<f_t> gr_refs(
          d_target_d_site_frac_.ref(),
          d_target_d_u_iso_.ref(),
          d_target_d_u_star_.ref(),
          d_target_d_occupancy_.ref(),
          d_target_d_fp_.ref(),
          d_target_d_fdp_.ref());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          gradients_direct_one_h<CosSinType, ScattererType>(
            cos_sin,
            unit_cell,
            space_group,
            miller_indices[i],
            scatterers,
            scattering_dict,
            d_target_d_f_calc[i],
            grad_flags,
            gr_refs);
        }
        {
          af::shared<std::size_t>
            perm = scattering_dict.scatterer_permutation();
          if (grad_flags.site) {
            d_target_d_site_frac_ = detail::unscramble(
              gr_refs.site, perm.const_ref(),
              static_cast<f_t>(scitbx::constants::two_pi));
          }
          if (grad_flags.u_iso) {
            d_target_d_u_iso_ = detail::unscramble(
              gr_refs.u_iso, perm.const_ref(),
              static_cast<f_t>(-scitbx::constants::two_pi_sq));
          }
          if (grad_flags.u_aniso) {
            d_target_d_u_star_ = detail::unscramble(
              gr_refs.u_star, perm.const_ref(),
              static_cast<f_t>(-scitbx::constants::two_pi_sq));
          }
          if (grad_flags.occupancy) {
            d_target_d_occupancy_ = detail::unscramble(
              gr_refs.occupancy, perm.const_ref());
          }
          if (grad_flags.fp) {
            d_target_d_fp_ = detail::unscramble(
              gr_refs.fp, perm.const_ref());
          }
          if (grad_flags.fdp) {
            d_target_d_fdp_ = detail::unscramble(
              gr_refs.fdp, perm.const_ref());
          }
        }
        if (n_parameters != 0) {
          gr_refs = detail::gradient_refs<f_t>(
            d_target_d_site_frac_.ref(),
            d_target_d_u_iso_.ref(),
            d_target_d_u_star_.ref(),
            d_target_d_occupancy_.ref(),
            d_target_d_fp_.ref(),
            d_target_d_fdp_.ref());
          BOOST_STATIC_ASSERT(packing_order_convention == 1);
          for(std::size_t i=0;i<scatterers.size();i++) {
            ScattererType const& scatterer = scatterers[i];
            if (grad_flags.site) {
              scitbx::vec3<f_t> d_target_d_site_cart =
                gr_refs.site[i] * unit_cell.fractionalization_matrix();
              for(std::size_t j=0;j<3;j++) {
                packed_.push_back(d_target_d_site_cart[j]);
              }
            }
            if (!scatterer.anisotropic_flag) {
              if (grad_flags.u_iso) {
                packed_.push_back(gr_refs.u_iso[i]);
              }
            }
            else {
              if (grad_flags.u_aniso) {
                scitbx::sym_mat3<double> d_target_d_u_cart =
                  adptbx::grad_u_star_as_u_cart(
                    unit_cell, gr_refs.u_star[i]);
                for(std::size_t j=0;j<6;j++) {
                  packed_.push_back(d_target_d_u_cart[j]);
                }
              }
            }
            if (grad_flags.occupancy) {
              packed_.push_back(gr_refs.occupancy[i]);
            }
            if (grad_flags.fp) {
              packed_.push_back(gr_refs.fp[i]);
            }
            if (grad_flags.fdp) {
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
