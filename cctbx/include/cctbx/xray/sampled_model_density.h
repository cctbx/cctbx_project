/* Based on phaser/src/MapFFT.cc by Airlie McCoy */

#ifndef CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H
#define CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H

#include <cctbx/xray/sampling_base.h>
#include <cctbx/maptbx/grid_tags.h>

namespace cctbx { namespace xray {

  template <typename FloatType=double,
            typename XrayScattererType=scatterer<> >
  class sampled_model_density
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

      sampled_model_density() {}

      sampled_model_density(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        grid_point_type const& fft_n_real,
        grid_point_type const& fft_m_real,
        FloatType const& u_extra=0.25,
        FloatType const& wing_cutoff=1.e-3,
        FloatType const& exp_table_one_over_step_size=-100,
        bool force_complex=false,
        bool electron_density_must_be_positive=true);

      real_map_type
      real_map() { return real_map_; }

      complex_map_type
      complex_map() { return complex_map_; }

      template <typename TagType>
      void
      apply_symmetry(maptbx::grid_tags<TagType> const& tags)
      {
        if (real_map_.size()) {
          tags.sum_sym_equiv_points(real_map_.ref());
        }
        else {
          tags.sum_sym_equiv_points(complex_map_.ref());
        }
      }

      void
      eliminate_u_extra_and_normalize(
        af::const_ref<miller::index<> > const& miller_indices,
        af::ref<std::complex<FloatType> > const& structure_factors) const
      {
        FloatType norm = this->unit_cell_.volume()
                       / this->map_accessor_.focus_size_1d();
        eliminate_u_extra(
          this->unit_cell_, this->u_extra_,
          miller_indices, structure_factors, norm);
      }

    private:
      real_map_type real_map_;
      complex_map_type complex_map_;
  };

  template <typename FloatType,
            typename XrayScattererType>
  sampled_model_density<FloatType, XrayScattererType>
  ::sampled_model_density(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<XrayScattererType> const& scatterers,
    grid_point_type const& fft_n_real,
    grid_point_type const& fft_m_real,
    FloatType const& u_extra,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    bool force_complex,
    bool electron_density_must_be_positive)
  :
    base_t(unit_cell, scatterers, u_extra, wing_cutoff,
           exp_table_one_over_step_size)
  {
    FloatType* map_begin;
    if (this->n_anomalous_scatterers_ == 0 && !force_complex) {
      this->map_accessor_ = accessor_type(fft_m_real, fft_n_real);
      real_map_.resize(this->map_accessor_);
      map_begin = real_map_.begin();
    }
    else {
      this->anomalous_flag_ = true;
      this->map_accessor_ = accessor_type(fft_n_real, fft_n_real);
      complex_map_.resize(this->map_accessor_);
      map_begin = reinterpret_cast<FloatType*>(complex_map_.begin());
    }
    grid_point_type const& grid_f = this->map_accessor_.focus();
    grid_point_type const& grid_a = this->map_accessor_.all();
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
      detail::caasf_fourier_transformed<FloatType, caasf_type> caasf_ft(
        exp_table,
        scatterer->caasf, scatterer->fp_fdp, scatterer->weight(),
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
        caasf_ft = detail::caasf_fourier_transformed<FloatType,caasf_type>(
          exp_table,
          scatterer->caasf, scatterer->fp_fdp, scatterer->weight(),
          u_cart, this->u_extra_);
      }
      grid_point_type pivot = detail::calc_nearest_grid_point(
        coor_frac, grid_f);
      // highly hand-optimized loop over points in shell
      grid_point_type g_min = pivot - shell.radii;
      grid_point_type g_max = pivot + shell.radii;
      grid_point_type gp;
      grid_point_element_type g01, g0112;
      FloatType f0, f1, f2;
      FloatType c00, c01, c11;
      for(gp[0] = g_min[0]; gp[0] <= g_max[0]; gp[0]++) {
        g01 = math::mod_positive(gp[0],grid_f[0]) * grid_a[1];
        f0 = FloatType(gp[0]) / grid_f[0] - coor_frac[0];
        c00 = orth_mx[0] * f0;
      for(gp[1] = g_min[1]; gp[1] <= g_max[1]; gp[1]++) {
        g0112 = (g01+math::mod_positive(gp[1],grid_f[1])) * grid_a[2];
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
        std::size_t i_map = g0112 + math::mod_positive(gp[2],grid_f[2]);
        if (this->anomalous_flag_) i_map *= 2;
        if (!scatterer->anisotropic_flag) {
          map_begin[i_map] += caasf_ft.rho_real(d_sq);
          if (fdp) {
            map_begin[i_map+1] += caasf_ft.rho_imag(d_sq);
          }
        }
        else {
          map_begin[i_map] += caasf_ft.rho_real(d);
          if (fdp) {
            map_begin[i_map+1] += caasf_ft.rho_imag(d);
          }
        }
      }}}
    }
    this->exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H
