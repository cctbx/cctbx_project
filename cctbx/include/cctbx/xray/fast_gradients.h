#ifndef CCTBX_XRAY_FAST_GRADIENTS_H
#define CCTBX_XRAY_FAST_GRADIENTS_H

#include <cctbx/xray/sampling_base.h>
#include <cctbx/xray/gradient_flags.h>

namespace cctbx { namespace xray {

  namespace detail {

    template <typename FloatType,
              typename CaasfType>
    class d_caasf_fourier_transformed
    :
      public caasf_fourier_transformed<FloatType, CaasfType>
    {
      public:
        d_caasf_fourier_transformed() {}

        d_caasf_fourier_transformed(
          gradient_flags const& gf,
          exponent_table<FloatType>& exp_table,
          CaasfType const& caasf,
          std::complex<FloatType> const& fp_fdp,
          FloatType const& occupancy,
          FloatType const& weight_without_occupancy,
          FloatType const& w,
          FloatType const& u_iso,
          FloatType const& u_extra)
        :
          caasf_fourier_transformed<FloatType, CaasfType>(
            exp_table, caasf, fp_fdp, w, u_iso, u_extra)
        {
          if (gf.u_iso || gf.occupancy || gf.fp || gf.fdp) {
            FloatType b_incl_extra = adptbx::u_as_b(u_iso + u_extra);
            if (gf.u_iso)
            {
              std::size_t i = 0;
              for(;i<caasf.n_ab();i++) {
                b_[i] = caasf.b(i) + b_incl_extra;
              }
              b_[i] = b_incl_extra;
            }
            if (gf.occupancy) {
              std::size_t i = 0;
              for(;i<caasf.n_ab();i++) {
                as_occupancy_real_[i]=isotropic_3d_gaussian_fourier_transform(
                  weight_without_occupancy * caasf.a(i),
                  caasf.b(i) + b_incl_extra);
              }
              as_occupancy_real_[i]=isotropic_3d_gaussian_fourier_transform(
                weight_without_occupancy * (caasf.c() + fp_fdp.real()),
                b_incl_extra);
              as_occupancy_imag_=isotropic_3d_gaussian_fourier_transform(
                weight_without_occupancy * fp_fdp.imag(),
                b_incl_extra);
            }
            if (gf.fp || gf.fdp) {
              FloatType d = b_incl_extra * b_incl_extra * b_incl_extra;
              eight_pi_pow_3_2_w_d_ = eight_pi_pow_3_2 * w / std::sqrt(d);
            }
          }
        }

        d_caasf_fourier_transformed(
          gradient_flags const& gf,
          exponent_table<FloatType>& exp_table,
          CaasfType const& caasf,
          std::complex<FloatType> const& fp_fdp,
          FloatType const& occupancy,
          FloatType const& weight_without_occupancy,
          FloatType const& w,
          scitbx::sym_mat3<FloatType> const& u_cart,
          FloatType const& u_extra)
        :
          caasf_fourier_transformed<FloatType, CaasfType>(
            exp_table, caasf, fp_fdp, w, u_cart, u_extra)
        {
          if (gf.u_aniso) {
            for(std::size_t i=0;i<caasf.n_ab();i++) {
              scitbx::sym_mat3<FloatType>
                b_all = compose_anisotropic_b_all(caasf.b(i), u_extra, u_cart);
              detb_[i] = b_all.determinant();
              bcfmt_[i] = b_all.co_factor_matrix_transposed();
            }
          }
          if (gf.u_aniso || gf.fp || gf.fdp) {
            scitbx::sym_mat3<FloatType>
              b_all = compose_anisotropic_b_all(0, u_extra, u_cart);
            FloatType d = b_all.determinant();
            if (gf.u_aniso) {
              static const std::size_t i = this->n_rho_real_terms - 1;
              detb_[i] = d;
              bcfmt_[i] = b_all.co_factor_matrix_transposed();
            }
            if (gf.fp || gf.fdp) {
              eight_pi_pow_3_2_w_d_ = eight_pi_pow_3_2 * w / std::sqrt(d);
            }
          }
          if (gf.occupancy) {
            std::size_t i = 0;
            for(;i<caasf.n_ab();i++) {
              as_occupancy_real_[i] =
                anisotropic_3d_gaussian_fourier_transform(
                  weight_without_occupancy * caasf.a(i),
                  compose_anisotropic_b_all(caasf.b(i), u_extra, u_cart));
            }
            as_occupancy_real_[i] =
              anisotropic_3d_gaussian_fourier_transform(
                weight_without_occupancy * (caasf.c() + fp_fdp.real()),
                compose_anisotropic_b_all(0, u_extra, u_cart));
            as_occupancy_imag_ =
              anisotropic_3d_gaussian_fourier_transform(
                weight_without_occupancy * (fp_fdp.imag()),
                compose_anisotropic_b_all(0, u_extra, u_cart));
          }
        }

        /* Mathematica script used to determine analytical gradients:
             r = as Exp[bs (d0^2+d1^2+d2^2)]
             g = {D[r,d0],D[r,d1],D[r,d2]}/r
           Conclusion:
             {g0,g1,g2} = 2 bs {d0,d1,d2} r
         */
        scitbx::vec3<FloatType>
        d_rho_real_d_site(scitbx::vec3<FloatType> const& d,
                          FloatType const& d_sq) const
        {
          scitbx::vec3<FloatType> drds(0,0,0);
          for (std::size_t i=0;i<this->n_rho_real_terms;i++) {
            drds += d*(this->bs_real_[i]*2*this->rho_real_term(d_sq, i));
          }
          return drds;
        }

        scitbx::vec3<FloatType>
        d_rho_imag_d_site(scitbx::vec3<FloatType> const& d,
                          FloatType const& d_sq) const
        {
          return d*(this->bs_imag_*2*this->rho_imag(d_sq));
        }

        scitbx::vec3<FloatType>
        d_rho_real_d_site(scitbx::vec3<FloatType> const& d) const
        {
          scitbx::vec3<FloatType> drds(0,0,0);
          for (std::size_t i=0;i<this->n_rho_real_terms;i++) {
            drds += this->aniso_bs_real_[i]*d*(2*this->rho_real_term(d, i));
          }
          return drds;
        }

        scitbx::vec3<FloatType>
        d_rho_imag_d_site(scitbx::vec3<FloatType> const& d) const
        {
          return this->aniso_bs_imag_*d*(2*this->rho_imag(d));
        }

        /* Mathematica script used to determine analytical gradients:
             f=a*(4 Pi / b)^(3/2) Exp[-4 Pi^2 d^2 / b]
             as=a*(4 Pi / b)^(3/2)
             bs=-4 Pi^2 / b
             FullSimplify[D[f,b]/f]
         */
        static
        FloatType
        d_rho_d_b_iso_term(
          FloatType const& d_sq,
          FloatType const& rho_term,
          FloatType const& b)
        {
          static const FloatType eight_pi_sq = 8 * scitbx::constants::pi_sq;
          return (3*b - eight_pi_sq*d_sq) / (2*b*b) * rho_term;
        }

        FloatType
        d_rho_real_d_b_iso(FloatType const& d_sq) const
        {
          FloatType drdb(0);
          for(std::size_t i=0;i<b_.size();i++) {
            drdb += d_rho_d_b_iso_term(
              d_sq, this->rho_real_term(d_sq, i), b_[i]);
          }
          return drdb;
        }

        FloatType
        d_rho_imag_d_b_iso(FloatType const& d_sq) const
        {
          return d_rho_d_b_iso_term(
            d_sq, this->rho_imag(d_sq), b_[this->n_rho_real_terms-1]);
        }

        /* Mathematica script used to determine analytical gradients:
             bcart={{b00+bx,b01,b02},{b01,b11+bx,b12},{b02,b12,b22+bx}}
             binv=Inverse[bcart]
             detb=Det[bcart]
             d={dx,dy,dz}
             r=ca/Sqrt[detb] Exp[cb d.binv.d]
             g={{D[r,b00], D[r,b01], D[r,b02]},
                {D[r,b01], D[r,b11], D[r,b12]},
                {D[r,b02], D[r,b12], D[r,b22]}}
             bcfmt=binv*detb
             bd=bcfmt.d
             cbd=cb/detb
             chk={{  cbd bd[[1]] bd[[1]] + bcfmt[[1,1]]/2,
                   2 cbd bd[[1]] bd[[2]] + bcfmt[[1,2]],
                   2 cbd bd[[1]] bd[[3]] + bcfmt[[1,3]]},
                  {2 cbd bd[[1]] bd[[2]] + bcfmt[[1,2]],
                     cbd bd[[2]] bd[[2]] + bcfmt[[2,2]]/2,
                   2 cbd bd[[2]] bd[[3]] + bcfmt[[2,3]]},
                  {2 cbd bd[[1]] bd[[3]] + bcfmt[[1,3]],
                   2 cbd bd[[2]] bd[[3]] + bcfmt[[2,3]],
                     cbd bd[[3]] bd[[3]] + bcfmt[[3,3]]/2}}/(-detb)
             FullSimplify[g-r*chk]
        */
        static
        scitbx::sym_mat3<FloatType>
        d_rho_d_b_cart_term(
          scitbx::vec3<FloatType> const& d,
          FloatType const& rho_term,
          FloatType const& detb,
          scitbx::sym_mat3<FloatType> const& bcfmt)
        {
          scitbx::sym_mat3<FloatType> drdb;
          scitbx::vec3<FloatType> bd = bcfmt * d;
          FloatType cbd = -four_pi_sq / detb;
          FloatType rd = rho_term / detb;
          for(std::size_t i=0;i<3;i++) {
            drdb[i] = rd * (cbd*bd[i]*bd[i] + bcfmt[i]*.5);
          }
          cbd *= 2;
          drdb[3] = rd * (cbd*bd[0]*bd[1] + bcfmt[3]);
          drdb[4] = rd * (cbd*bd[0]*bd[2] + bcfmt[4]);
          drdb[5] = rd * (cbd*bd[1]*bd[2] + bcfmt[5]);
          return drdb;
        }

        scitbx::sym_mat3<FloatType>
        d_rho_real_d_b_cart(scitbx::vec3<FloatType> const& d) const
        {
          scitbx::sym_mat3<FloatType> drdb(0,0,0,0,0,0);
          for (std::size_t i=0;i<this->n_rho_real_terms;i++) {
            drdb += d_rho_d_b_cart_term(
              d, this->rho_real_term(d, i), detb_[i], bcfmt_[i]);
          }
          return drdb;
        }

        scitbx::sym_mat3<FloatType>
        d_rho_imag_d_b_cart(scitbx::vec3<FloatType> const& d) const
        {
          return d_rho_d_b_cart_term(
            d, this->rho_imag(d),
            detb_[this->n_rho_real_terms-1],
            bcfmt_[this->n_rho_real_terms-1]);
        }

        template <typename DistanceType>
        FloatType
        d_rho_real_d_occupancy(DistanceType const& d_or_d_sq) const
        {
          FloatType drdo(0);
          for (std::size_t i=0;i<this->n_rho_real_terms;i++) {
            drdo -= as_occupancy_real_[i] * this->exp_term(d_or_d_sq, i);
          }
          return drdo;
        }

        template <typename DistanceType>
        FloatType
        d_rho_imag_d_occupancy(DistanceType const& d_or_d_sq) const
        {
          return -as_occupancy_imag_ * this->exp_term(d_or_d_sq);
        }

        template <typename DistanceType>
        FloatType
        d_rho_real_d_fp(DistanceType const& d_or_d_sq) const
        {
          return -eight_pi_pow_3_2_w_d_
                 * this->exp_term(d_or_d_sq, this->n_rho_real_terms-1);
        }

        template <typename DistanceType>
        FloatType
        d_rho_imag_d_fdp(DistanceType const& d_or_d_sq) const
        {
          return d_rho_real_d_fp(d_or_d_sq);
        }

      protected:
        af::tiny<FloatType, CaasfType::n_plus_1> b_;
        af::tiny<FloatType, CaasfType::n_plus_1> detb_;
        af::tiny<scitbx::sym_mat3<FloatType>, CaasfType::n_plus_1> bcfmt_;
        af::tiny<FloatType, CaasfType::n_plus_1> as_occupancy_real_;
        FloatType as_occupancy_imag_;
        FloatType eight_pi_pow_3_2_w_d_;
    };

  } // namespace detail

  template <typename FloatType=double,
            typename XrayScattererType=scatterer<> >
  class fast_gradients
  :
    public sampling_base<FloatType, XrayScattererType>
  {
    public:
      typedef sampling_base<FloatType, XrayScattererType> base_t;
      typedef typename base_t::xray_scatterer_type xray_scatterer_type;
      typedef typename base_t::caasf_type caasf_type;
      typedef typename base_t::accessor_type accessor_type;
      typedef typename base_t::grid_point_type grid_point_type;
      typedef typename base_t::grid_point_element_type grid_point_element_type;
      typedef typename base_t::real_map_type real_map_type;
      typedef typename base_t::complex_map_type complex_map_type;

      fast_gradients() {}

      fast_gradients(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        af::const_ref<std::complex<FloatType>, accessor_type> const&
          ft_d_target_d_f_calc,
        gradient_flags const& grad_flags,
        FloatType const& u_extra=0.25,
        FloatType const& wing_cutoff=1.e-3,
        FloatType const& exp_table_one_over_step_size=-100,
        bool electron_density_must_be_positive=true);

      af::shared<scitbx::vec3<FloatType> >
      d_target_d_site_cart() const { return d_target_d_site_cart_; }

      af::shared<FloatType>
      d_target_d_u_iso() const { return d_target_d_u_iso_; }

      af::shared<scitbx::sym_mat3<FloatType> >
      d_target_d_u_cart() const { return d_target_d_u_cart_; }

      af::shared<FloatType>
      d_target_d_occupancy() const { return d_target_d_occupancy_; }

      af::shared<FloatType>
      d_target_d_fp() const { return d_target_d_fp_; }

      af::shared<FloatType>
      d_target_d_fdp() const { return d_target_d_fdp_; }

    private:
      af::shared<scitbx::vec3<FloatType> > d_target_d_site_cart_;
      af::shared<FloatType> d_target_d_u_iso_;
      af::shared<scitbx::sym_mat3<FloatType> > d_target_d_u_cart_;
      af::shared<FloatType> d_target_d_occupancy_;
      af::shared<FloatType> d_target_d_fp_;
      af::shared<FloatType> d_target_d_fdp_;
  };

  template <typename FloatType,
            typename XrayScattererType>
  fast_gradients<FloatType, XrayScattererType>
  ::fast_gradients(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<XrayScattererType> const& scatterers,
    af::const_ref<std::complex<FloatType>, accessor_type> const&
      ft_d_target_d_f_calc,
    gradient_flags const& grad_flags,
    FloatType const& u_extra,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    bool electron_density_must_be_positive)
  :
    base_t(unit_cell, scatterers, u_extra, wing_cutoff,
           exp_table_one_over_step_size)
  {
    if (this->n_anomalous_scatterers_ != 0) {
      this->anomalous_flag_ = true;
    }
    this->map_accessor_ = ft_d_target_d_f_calc.accessor();
    grid_point_type const& grid_f = this->map_accessor_.focus();
    detail::exponent_table<FloatType> exp_table(exp_table_one_over_step_size);
    scitbx::mat3<FloatType>
      orth_mx = this->unit_cell_.orthogonalization_matrix();
    const xray_scatterer_type* scatterer;
    for(scatterer=scatterers.begin();scatterer!=scatterers.end();scatterer++)
    {
      if (scatterer->weight() == 0) continue;
      FloatType fdp = scatterer->fp_fdp.imag();
      fractional<FloatType> coor_frac = scatterer->site;
      FloatType u_iso;
      scitbx::sym_mat3<FloatType> u_cart;
      if (!scatterer->anisotropic_flag) {
        u_iso = scatterer->u_iso;
      }
      else {
        u_cart = adptbx::u_star_as_u_cart(this->unit_cell_, scatterer->u_star);
        scitbx::vec3<FloatType> ev = adptbx::eigenvalues(u_cart);
        CCTBX_ASSERT(adptbx::is_positive_definite(ev));
        u_iso = af::max(ev);
      }
      CCTBX_ASSERT(u_iso >= 0);
      detail::d_caasf_fourier_transformed<FloatType, caasf_type> caasf_ft(
        grad_flags,
        exp_table,
        scatterer->caasf, scatterer->fp_fdp, scatterer->occupancy,
        scatterer->weight_without_occupancy(), scatterer->weight(),
        u_iso, this->u_extra_);
      detail::calc_shell<FloatType, grid_point_type> shell(
        this->unit_cell_, this->wing_cutoff_, grid_f, caasf_ft);
      detail::array_update_max(this->max_shell_radii_, shell.radii);
      if (electron_density_must_be_positive) {
        if (   caasf_ft.rho_real_0() < 0
            || caasf_ft.rho_real(shell.max_d_sq) < 0) {

          throw error("Negative electron density at sampling point.");
        }
      }
      if (scatterer->anisotropic_flag) {
        caasf_ft = detail::d_caasf_fourier_transformed<FloatType,caasf_type>(
          grad_flags,
          exp_table,
          scatterer->caasf, scatterer->fp_fdp, scatterer->occupancy,
          scatterer->weight_without_occupancy(), scatterer->weight(),
          u_cart, this->u_extra_);
      }
      scitbx::vec3<FloatType> gr_site(0,0,0);
      FloatType gr_b_iso(0);
      scitbx::sym_mat3<FloatType> gr_b_cart(0,0,0,0,0,0);
      FloatType gr_occupancy(0);
      FloatType gr_fp(0);
      FloatType gr_fdp(0);
      grid_point_type pivot = detail::calc_nearest_grid_point(
        coor_frac, grid_f);
      // highly hand-optimized loop over points in shell
      grid_point_type g_min = pivot - shell.radii;
      grid_point_type g_max = pivot + shell.radii;
      grid_point_type gp;
      FloatType f0, f1, f2;
      FloatType c00, c01, c11;
      for(gp[0] = g_min[0]; gp[0] <= g_max[0]; gp[0]++) {
        f0 = FloatType(gp[0]) / grid_f[0] - coor_frac[0];
        c00 = orth_mx[0] * f0;
      for(gp[1] = g_min[1]; gp[1] <= g_max[1]; gp[1]++) {
        f1 = FloatType(gp[1]) / grid_f[1] - coor_frac[1];
        c01 = orth_mx[1] * f1 + c00;
        c11 = orth_mx[4] * f1;
      for(gp[2] = g_min[2]; gp[2] <= g_max[2]; gp[2]++) {
        f2 = FloatType(gp[2]) / grid_f[2] - coor_frac[2];
        scitbx::vec3<FloatType> d(
          orth_mx[2] * f2 + c01,
          orth_mx[5] * f2 + c11,
          orth_mx[8] * f2);
        FloatType d_sq = d.length_sq();
        if (d_sq > shell.max_d_sq) continue;
        std::complex<FloatType> ft_dt_dfc_gp = ft_d_target_d_f_calc(gp);
        FloatType f_real = -ft_dt_dfc_gp.real();
        FloatType f_imag = -ft_dt_dfc_gp.imag();
        if (grad_flags.site) {
          if (!scatterer->anisotropic_flag) {
            gr_site += f_real * caasf_ft.d_rho_real_d_site(d, d_sq);
            if (fdp) {
              gr_site += f_imag * caasf_ft.d_rho_imag_d_site(d, d_sq);
            }
          }
          else {
            gr_site += f_real * caasf_ft.d_rho_real_d_site(d);
            if (fdp) {
              gr_site += f_imag * caasf_ft.d_rho_imag_d_site(d);
            }
          }
        }
        if (!scatterer->anisotropic_flag) {
          if (grad_flags.u_iso) {
            gr_b_iso += f_real * caasf_ft.d_rho_real_d_b_iso(d_sq);
            if (fdp) {
              gr_b_iso += f_imag * caasf_ft.d_rho_imag_d_b_iso(d_sq);
            }
          }
        }
        else {
          if (grad_flags.u_aniso) {
            gr_b_cart += f_real * caasf_ft.d_rho_real_d_b_cart(d);
            if (fdp) {
              gr_b_cart += f_imag * caasf_ft.d_rho_imag_d_b_cart(d);
            }
          }
        }
        if (grad_flags.occupancy) {
          if (!scatterer->anisotropic_flag) {
            gr_occupancy += f_real * caasf_ft.d_rho_real_d_occupancy(d_sq);
            if (fdp) {
              gr_occupancy += f_imag * caasf_ft.d_rho_imag_d_occupancy(d_sq);
            }
          }
          else {
            gr_occupancy += f_real * caasf_ft.d_rho_real_d_occupancy(d);
            if (fdp) {
              gr_occupancy += f_imag * caasf_ft.d_rho_imag_d_occupancy(d);
            }
          }
        }
        if (grad_flags.fp) {
          if (!scatterer->anisotropic_flag) {
            gr_fp += f_real * caasf_ft.d_rho_real_d_fp(d_sq);
          }
          else {
            gr_fp += f_real * caasf_ft.d_rho_real_d_fp(d);
          }
        }
        if (grad_flags.fdp) {
          if (!scatterer->anisotropic_flag) {
            gr_fdp += f_imag * caasf_ft.d_rho_imag_d_fdp(d_sq);
          }
          else {
            gr_fdp += f_imag * caasf_ft.d_rho_imag_d_fdp(d);
          }
        }
      }}}
      if (grad_flags.site) {
        d_target_d_site_cart_.push_back(gr_site);
      }
      if (grad_flags.u_iso) {
        d_target_d_u_iso_.push_back(adptbx::u_as_b(gr_b_iso));
      }
      if (grad_flags.u_aniso) {
        d_target_d_u_cart_.push_back(adptbx::u_as_b(gr_b_cart));
      }
      if (grad_flags.occupancy) {
        d_target_d_occupancy_.push_back(gr_occupancy);
      }
      if (grad_flags.fp) {
        d_target_d_fp_.push_back(gr_fp);
      }
      if (grad_flags.fdp) {
        d_target_d_fdp_.push_back(gr_fdp);
      }
    }
    this->exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_FAST_GRADIENTS_H
