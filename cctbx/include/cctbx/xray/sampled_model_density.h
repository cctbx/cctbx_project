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
      typedef typename base_t::accessor_type accessor_type;
      typedef typename base_t::grid_point_type grid_point_type;
      typedef typename base_t::grid_point_element_type grid_point_element_type;
      typedef typename base_t::real_map_type real_map_type;
      typedef typename base_t::complex_map_type complex_map_type;

      sampled_model_density() {}

      sampled_model_density(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        scattering_dictionary const& scattering_dict,
        grid_point_type const& fft_n_real,
        grid_point_type const& fft_m_real,
        FloatType const& u_extra=0.25,
        FloatType const& wing_cutoff=1.e-6,
        FloatType const& exp_table_one_over_step_size=-100,
        bool force_complex=false,
        bool electron_density_must_be_positive=true,
        FloatType const& tolerance_positive_definite=1.e-5);

      real_map_type
      real_map() { return real_map_; }

      complex_map_type
      complex_map() { return complex_map_; }

      void
      eliminate_u_extra_and_normalize(
        af::const_ref<miller::index<> > const& miller_indices,
        af::ref<std::complex<FloatType> > const& structure_factors) const
      {
        FloatType multiplier = this->unit_cell_.volume()
                             / this->map_accessor_.focus_size_1d();
        apply_u_extra(
          this->unit_cell_, this->u_extra_,
          miller_indices, structure_factors, multiplier);
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
    scattering_dictionary const& scattering_dict,
    grid_point_type const& fft_n_real,
    grid_point_type const& fft_m_real,
    FloatType const& u_extra,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    bool force_complex,
    bool electron_density_must_be_positive,
    FloatType const& tolerance_positive_definite)
  :
    base_t(unit_cell, scatterers, u_extra, wing_cutoff,
           exp_table_one_over_step_size, tolerance_positive_definite)
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
    typedef scattering_dictionary::dict_type dict_type;
    typedef dict_type::const_iterator dict_iter;
    dict_type const& scd = scattering_dict.dict();
    for(dict_iter di=scd.begin();di!=scd.end();di++) {
      eltbx::caasf::custom const& caasf = di->second.coefficients;
      af::const_ref<std::size_t>
        member_indices = di->second.member_indices.const_ref();
      for(std::size_t mi=0;mi<member_indices.size();mi++) {
        XrayScattererType const& scatterer = scatterers[member_indices[mi]];
        CCTBX_ASSERT(scatterer.weight() >= 0);
        if (scatterer.weight() == 0) continue;
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
        detail::caasf_fourier_transformed<FloatType> caasf_ft(
          exp_table,
          caasf, scatterer.fp, scatterer.fdp, scatterer.weight(),
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
        if (scatterer.anisotropic_flag) {
          caasf_ft = detail::caasf_fourier_transformed<FloatType>(
            exp_table,
            caasf, scatterer.fp, scatterer.fdp, scatterer.weight(),
            u_cart, this->u_extra_);
        }
        grid_point_type pivot = detail::calc_nearest_grid_point(
          coor_frac, grid_f);
#       include <cctbx/xray/sampling_loop.h>
          if (this->anomalous_flag_) i_map *= 2;
          if (!scatterer.anisotropic_flag) {
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
    }
    this->exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H
