#ifndef CCTBX_XRAY_FAST_GRADIENTS_H
#define CCTBX_XRAY_FAST_GRADIENTS_H

#include <cctbx/xray/sampling_base.h>
#include <cctbx/xray/gradient_flags.h>
#include <cctbx/xray/packing_order.h>

namespace cctbx { namespace xray {

  namespace detail {

    template <typename FloatType>
    class d_gaussian_fourier_transformed
    :
      public gaussian_fourier_transformed<FloatType>
    {
      public:
        typedef gaussian_fourier_transformed<FloatType> base_t;

        d_gaussian_fourier_transformed() {}

        d_gaussian_fourier_transformed(
          gradient_flags const& gf,
          exponent_table<FloatType>& exp_table,
          eltbx::xray_scattering::gaussian const& gaussian,
          FloatType const& fp,
          FloatType const& fdp,
          FloatType const& occupancy,
          FloatType const& weight_without_occupancy,
          FloatType const& w,
          FloatType const& u_iso,
          FloatType const& u_extra)
        :
          base_t(exp_table, gaussian, fp, fdp, w, u_iso, u_extra),
          i_const_term(gaussian.n_ab())
        {
          if (gf.u_iso || gf.occupancy || gf.fp || gf.fdp) {
            FloatType b_incl_extra = adptbx::u_as_b(u_iso + u_extra);
            if (gf.u_iso)
            {
              std::size_t i = 0;
              for(;i<gaussian.n_ab();i++) {
                b_[i] = gaussian.b(i) + b_incl_extra;
              }
              b_[i] = b_incl_extra;
            }
            if (gf.occupancy) {
              std::size_t i = 0;
              for(;i<gaussian.n_ab();i++) {
                as_occupancy_real_[i]=isotropic_3d_gaussian_fourier_transform(
                  weight_without_occupancy * gaussian.a(i),
                  gaussian.b(i) + b_incl_extra);
              }
              if (this->n_rho_real_terms > gaussian.n_ab()) {
                as_occupancy_real_[i]=isotropic_3d_gaussian_fourier_transform(
                  weight_without_occupancy * (gaussian.c() + fp),
                  b_incl_extra);
              }
              as_occupancy_imag_=isotropic_3d_gaussian_fourier_transform(
                weight_without_occupancy * fdp,
                b_incl_extra);
            }
            if (gf.fp || gf.fdp) {
              FloatType d = b_incl_extra * b_incl_extra * b_incl_extra;
              eight_pi_pow_3_2_w_d_ = eight_pi_pow_3_2 * w / std::sqrt(d);
            }
          }
        }

        d_gaussian_fourier_transformed(
          gradient_flags const& gf,
          exponent_table<FloatType>& exp_table,
          eltbx::xray_scattering::gaussian const& gaussian,
          FloatType const& fp,
          FloatType const& fdp,
          FloatType const& occupancy,
          FloatType const& weight_without_occupancy,
          FloatType const& w,
          scitbx::sym_mat3<FloatType> const& u_cart,
          FloatType const& u_extra)
        :
          base_t(exp_table, gaussian, fp, fdp, w, u_cart, u_extra),
          i_const_term(gaussian.n_ab())
        {
          if (gf.u_aniso) {
            for(std::size_t i=0;i<gaussian.n_ab();i++) {
              scitbx::sym_mat3<FloatType>
                b_all = compose_anisotropic_b_all(
                  gaussian.b(i), u_extra, u_cart);
              detb_[i] = b_all.determinant();
              bcfmt_[i] = b_all.co_factor_matrix_transposed();
            }
          }
          if (gf.u_aniso || gf.fp || gf.fdp) {
            scitbx::sym_mat3<FloatType>
              b_all = compose_anisotropic_b_all(0, u_extra, u_cart);
            FloatType d = b_all.determinant();
            if (gf.u_aniso) {
              detb_[i_const_term] = d;
              bcfmt_[i_const_term] = b_all.co_factor_matrix_transposed();
            }
            if (gf.fp || gf.fdp) {
              eight_pi_pow_3_2_w_d_ = eight_pi_pow_3_2 * w / std::sqrt(d);
            }
          }
          if (gf.occupancy) {
            std::size_t i = 0;
            for(;i<gaussian.n_ab();i++) {
              as_occupancy_real_[i] =
                anisotropic_3d_gaussian_fourier_transform(
                  weight_without_occupancy * gaussian.a(i),
                  compose_anisotropic_b_all(gaussian.b(i), u_extra, u_cart));
            }
            if (this->n_rho_real_terms > gaussian.n_ab()) {
              as_occupancy_real_[i] =
                anisotropic_3d_gaussian_fourier_transform(
                  weight_without_occupancy * (gaussian.c() + fp),
                  compose_anisotropic_b_all(0, u_extra, u_cart));
            }
            as_occupancy_imag_ =
              anisotropic_3d_gaussian_fourier_transform(
                weight_without_occupancy * fdp,
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
          for(std::size_t i=0;i<this->n_rho_real_terms;i++) {
            drdb += d_rho_d_b_iso_term(
              d_sq, this->rho_real_term(d_sq, i), b_[i]);
          }
          return drdb;
        }

        FloatType
        d_rho_imag_d_b_iso(FloatType const& d_sq) const
        {
          return d_rho_d_b_iso_term(
            d_sq, this->rho_imag(d_sq), b_[i_const_term]);
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
            detb_[i_const_term],
            bcfmt_[i_const_term]);
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
                 * this->exp_term(d_or_d_sq, i_const_term);
        }

        template <typename DistanceType>
        FloatType
        d_rho_imag_d_fdp(DistanceType const& d_or_d_sq) const
        {
          return d_rho_real_d_fp(d_or_d_sq);
        }

      protected:
        std::size_t i_const_term;
        af::tiny<FloatType, base_t::max_n_rho_real_terms> b_;
        af::tiny<FloatType, base_t::max_n_rho_real_terms> detb_;
        af::tiny<scitbx::sym_mat3<FloatType>, base_t::max_n_rho_real_terms>
          bcfmt_;
        af::tiny<FloatType, base_t::max_n_rho_real_terms> as_occupancy_real_;
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
      typedef typename base_t::accessor_type accessor_type;
      typedef typename base_t::grid_point_type grid_point_type;
      typedef typename base_t::grid_point_element_type grid_point_element_type;
      typedef typename base_t::real_map_type real_map_type;
      typedef typename base_t::complex_map_type complex_map_type;

      fast_gradients() {}

      fast_gradients(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        scattering_dictionary const& scattering_dict,
        af::const_ref<FloatType, accessor_type> const&
          ft_d_target_d_f_calc_real,
        af::const_ref<std::complex<FloatType>, accessor_type> const&
          ft_d_target_d_f_calc_complex,
        gradient_flags const& grad_flags,
        std::size_t n_parameters=0,
        FloatType const& u_extra=0.25,
        FloatType const& wing_cutoff=1.e-6,
        FloatType const& exp_table_one_over_step_size=-100,
        bool electron_density_must_be_positive=true,
        FloatType const& tolerance_positive_definite=1.e-5);

      af::shared<FloatType>
      packed() const { return packed_; }

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
      af::shared<FloatType> packed_;
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
    scattering_dictionary const& scattering_dict,
    af::const_ref<FloatType, accessor_type> const&
      ft_d_target_d_f_calc_real,
    af::const_ref<std::complex<FloatType>, accessor_type> const&
      ft_d_target_d_f_calc_complex,
    gradient_flags const& grad_flags,
    std::size_t n_parameters,
    FloatType const& u_extra,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    bool electron_density_must_be_positive,
    FloatType const& tolerance_positive_definite)
  :
    base_t(unit_cell, scatterers, u_extra, wing_cutoff,
           exp_table_one_over_step_size, tolerance_positive_definite)
  {
    CCTBX_ASSERT(scattering_dict.n_scatterers() == scatterers.size());
    CCTBX_ASSERT(   (ft_d_target_d_f_calc_real.size() == 0)
                 != (ft_d_target_d_f_calc_complex.size() == 0));
    if (this->n_anomalous_scatterers_ != 0) {
      this->anomalous_flag_ = true;
    }
    if (n_parameters != 0) {
      packed_.reserve(n_parameters);
    }
    else {
      if (grad_flags.site) d_target_d_site_cart_.reserve(scatterers.size());
      if (grad_flags.u_iso) d_target_d_u_iso_.reserve(scatterers.size());
      if (grad_flags.u_aniso) d_target_d_u_cart_.reserve(scatterers.size());
      if (grad_flags.occupancy) d_target_d_occupancy_.reserve(
                                  scatterers.size());
      if (grad_flags.fp) d_target_d_fp_.reserve(scatterers.size());
      if (grad_flags.fdp) d_target_d_fdp_.reserve(scatterers.size());
    }
    bool gradient_map_is_real = ft_d_target_d_f_calc_real.size() != 0;
    if (gradient_map_is_real) {
      this->map_accessor_ = ft_d_target_d_f_calc_real.accessor();
    }
    else {
      this->map_accessor_ = ft_d_target_d_f_calc_complex.accessor();
    }
    grid_point_type const& grid_f = this->map_accessor_.focus();
    grid_point_type const& grid_a = this->map_accessor_.all();
    detail::exponent_table<FloatType> exp_table(exp_table_one_over_step_size);
    scitbx::mat3<FloatType>
      orth_mx = this->unit_cell_.orthogonalization_matrix();
    typedef scattering_dictionary::dict_type dict_type;
    typedef dict_type::const_iterator dict_iter;
    dict_type const& scd = scattering_dict.dict();
    for(dict_iter di=scd.begin();di!=scd.end();di++) {
      eltbx::xray_scattering::gaussian const& gaussian = di->second.gaussian;
      af::const_ref<std::size_t>
        member_indices = di->second.member_indices.const_ref();
      for(std::size_t mi=0;mi<member_indices.size();mi++) {
        XrayScattererType const& scatterer = scatterers[member_indices[mi]];
        scitbx::vec3<FloatType> gr_site(0,0,0);
        FloatType gr_b_iso(0);
        scitbx::sym_mat3<FloatType> gr_b_cart(0,0,0,0,0,0);
        FloatType gr_occupancy(0);
        FloatType gr_fp(0);
        FloatType gr_fdp(0);
        CCTBX_ASSERT(scatterer.weight() >= 0);
        if (scatterer.weight() != 0) {
          FloatType fdp = scatterer.fdp;
          fractional<FloatType> coor_frac = scatterer.site;
          FloatType u_iso;
          scitbx::sym_mat3<FloatType> u_cart;
          if (!scatterer.anisotropic_flag) {
            u_iso = scatterer.u_iso;
          }
          else {
            u_iso = this->get_u_cart_and_u_iso(scatterer.u_star, u_cart);
          }
          CCTBX_ASSERT(u_iso >= 0);
          detail::d_gaussian_fourier_transformed<FloatType> gaussian_ft(
            grad_flags,
            exp_table,
            gaussian, scatterer.fp, scatterer.fdp, scatterer.occupancy,
            scatterer.weight_without_occupancy(), scatterer.weight(),
            u_iso, this->u_extra_);
          detail::calc_shell<FloatType, grid_point_type> shell(
            this->unit_cell_, this->wing_cutoff_, grid_f, gaussian_ft);
          detail::array_update_max(this->max_shell_radii_, shell.radii);
          if (electron_density_must_be_positive) {
            if (   gaussian_ft.rho_real_0() < 0
                || gaussian_ft.rho_real(shell.max_d_sq) < 0) {

              throw error("Negative electron density at sampling point.");
            }
          }
          if (scatterer.anisotropic_flag) {
            gaussian_ft = detail::d_gaussian_fourier_transformed<FloatType>(
              grad_flags,
              exp_table,
              gaussian, scatterer.fp, scatterer.fdp, scatterer.occupancy,
              scatterer.weight_without_occupancy(), scatterer.weight(),
              u_cart, this->u_extra_);
          }
          grid_point_type pivot = detail::calc_nearest_grid_point(
            coor_frac, grid_f);
          FloatType f_real(-1.e20); // could be uninitialized
          FloatType f_imag(-1.e20); // could be uninitialized
#         include <cctbx/xray/sampling_loop.h>
            if (gradient_map_is_real) {
              f_real = -ft_d_target_d_f_calc_real[i_map];
            }
            else {
              std::complex<FloatType> const&
                ft_dt_dfc_gp = ft_d_target_d_f_calc_complex[i_map];
              f_real = -ft_dt_dfc_gp.real();
              f_imag = -ft_dt_dfc_gp.imag();
            }
            if (grad_flags.site) {
              if (!scatterer.anisotropic_flag) {
                gr_site += f_real * gaussian_ft.d_rho_real_d_site(d, d_sq);
                if (fdp && !gradient_map_is_real) {
                  gr_site += f_imag * gaussian_ft.d_rho_imag_d_site(d, d_sq);
                }
              }
              else {
                gr_site += f_real * gaussian_ft.d_rho_real_d_site(d);
                if (fdp && !gradient_map_is_real) {
                  gr_site += f_imag * gaussian_ft.d_rho_imag_d_site(d);
                }
              }
            }
            if (!scatterer.anisotropic_flag) {
              if (grad_flags.u_iso) {
                gr_b_iso += f_real * gaussian_ft.d_rho_real_d_b_iso(d_sq);
                if (fdp && !gradient_map_is_real) {
                  gr_b_iso += f_imag * gaussian_ft.d_rho_imag_d_b_iso(d_sq);
                }
              }
            }
            else {
              if (grad_flags.u_aniso) {
                gr_b_cart += f_real * gaussian_ft.d_rho_real_d_b_cart(d);
                if (fdp && !gradient_map_is_real) {
                  gr_b_cart += f_imag * gaussian_ft.d_rho_imag_d_b_cart(d);
                }
              }
            }
            if (grad_flags.occupancy) {
              if (!scatterer.anisotropic_flag) {
                gr_occupancy += f_real
                              * gaussian_ft.d_rho_real_d_occupancy(d_sq);
                if (fdp && !gradient_map_is_real) {
                  gr_occupancy += f_imag
                                * gaussian_ft.d_rho_imag_d_occupancy(d_sq);
                }
              }
              else {
                gr_occupancy += f_real
                              * gaussian_ft.d_rho_real_d_occupancy(d);
                if (fdp && !gradient_map_is_real) {
                  gr_occupancy += f_imag
                                * gaussian_ft.d_rho_imag_d_occupancy(d);
                }
              }
            }
            if (grad_flags.fp) {
              if (!scatterer.anisotropic_flag) {
                gr_fp += f_real * gaussian_ft.d_rho_real_d_fp(d_sq);
              }
              else {
                gr_fp += f_real * gaussian_ft.d_rho_real_d_fp(d);
              }
            }
            if (grad_flags.fdp && !gradient_map_is_real) {
              if (!scatterer.anisotropic_flag) {
                gr_fdp += f_imag * gaussian_ft.d_rho_imag_d_fdp(d_sq);
              }
              else {
                gr_fdp += f_imag * gaussian_ft.d_rho_imag_d_fdp(d);
              }
            }
          }}}
        }
        if (n_parameters != 0) {
          BOOST_STATIC_ASSERT(packing_order_convention == 1);
          if (grad_flags.site) {
            for(std::size_t i=0;i<3;i++) {
              packed_.push_back(gr_site[i]);
            }
          }
          if (!scatterer.anisotropic_flag) {
            if (grad_flags.u_iso) {
              packed_.push_back(adptbx::u_as_b(gr_b_iso));
            }
          }
          else {
            if (grad_flags.u_aniso) {
              for(std::size_t i=0;i<6;i++) {
                packed_.push_back(adptbx::u_as_b(gr_b_cart[i]));
              }
            }
          }
          if (grad_flags.occupancy) {
            packed_.push_back(gr_occupancy);
          }
          if (grad_flags.fp) {
            packed_.push_back(gr_fp);
          }
          if (grad_flags.fdp) {
            packed_.push_back(gr_fdp);
          }
        }
        else {
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
      }
    }
    if (n_parameters != 0) {
      CCTBX_ASSERT(packed_.size() == n_parameters);
    }
    this->exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_FAST_GRADIENTS_H
