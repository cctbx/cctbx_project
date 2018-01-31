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
      typedef typename base_t::u_radius_cache_t u_radius_cache_t;
      typedef typename base_t::u_cart_cache_t u_cart_cache_t;

      sampled_model_density() {}

      sampled_model_density(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        xray::scattering_type_registry const& scattering_type_registry,
        grid_point_type const& fft_n_real,
        grid_point_type const& fft_m_real,
        FloatType const& u_base=0.25,
        FloatType const& wing_cutoff=1.e-3,
        FloatType const& exp_table_one_over_step_size=-100,
        bool force_complex=false,
        bool sampled_density_must_be_positive=false,
        FloatType const& tolerance_positive_definite=1.e-5,
        bool use_u_base_as_u_extra=false,
        int store_grid_indices_for_each_scatterer=0);

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

      af::shared<std::vector<unsigned> > const&
      grid_indices_for_each_scatterer() const
      {
        return grid_indices_for_each_scatterer_;
      }

    protected:
      real_map_type real_map_;
      complex_map_type complex_map_;
      af::shared<std::vector<unsigned> > grid_indices_for_each_scatterer_;
  };

  template <typename FloatType,
            typename XrayScattererType>
  sampled_model_density<FloatType, XrayScattererType>
  ::sampled_model_density(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<XrayScattererType> const& scatterers,
    xray::scattering_type_registry const& scattering_type_registry,
    grid_point_type const& fft_n_real,
    grid_point_type const& fft_m_real,
    FloatType const& u_base,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    bool force_complex,
    bool sampled_density_must_be_positive,
    FloatType const& tolerance_positive_definite,
    bool use_u_base_as_u_extra,
    int store_grid_indices_for_each_scatterer)
  :
    base_t(unit_cell, scatterers, scattering_type_registry,
           u_base, wing_cutoff,
           exp_table_one_over_step_size, tolerance_positive_definite,
           use_u_base_as_u_extra)
  {
    FloatType* map_begin = 0;
    if (this->n_anomalous_scatterers_ == 0 && !force_complex) {
      this->map_accessor_ = accessor_type(fft_m_real, fft_n_real);
      if (store_grid_indices_for_each_scatterer >= 0) {
        real_map_.resize(this->map_accessor_);
        map_begin = real_map_.begin();
      }
    }
    else {
      this->anomalous_flag_ = true;
      this->map_accessor_ = accessor_type(fft_n_real, fft_n_real);
      if (store_grid_indices_for_each_scatterer >= 0) {
        complex_map_.resize(this->map_accessor_);
        map_begin = reinterpret_cast<FloatType*>(complex_map_.begin());
      }
    }
    std::vector<unsigned>* gifes = 0;
    if (store_grid_indices_for_each_scatterer != 0) {
      grid_indices_for_each_scatterer_.resize(scatterers.size());
      gifes = grid_indices_for_each_scatterer_.begin();
    }
    std::vector<unsigned> grid_indices_for_one_scatterer;
    grid_point_type const& grid_f = this->map_accessor_.focus();
    grid_point_type const& grid_a = this->map_accessor_.all();
    detail::exponent_table<FloatType> exp_table(exp_table_one_over_step_size);
    scitbx::mat3<FloatType>
      orth_mx = this->unit_cell_.orthogonalization_matrix();
    std::vector<int> gp1g;
    std::vector<FloatType> o1f1_;
    std::vector<FloatType> o4f1_;
    std::vector<int> gp2g;
    std::vector<FloatType> o2f2_;
    std::vector<FloatType> o5f2_;
    std::vector<FloatType> o8f2_;
    typename u_radius_cache_t::const_iterator u_radius =
      this->u_radius_cache_.begin();
    typename u_cart_cache_t::const_iterator u_cart =
      this->u_cart_cache_.begin();
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      XrayScattererType const& scatterer = scatterers[i_seq];
      if (scatterer.weight() == 0) continue;
      eltbx::xray_scattering::gaussian const&
        gaussian = scattering_type_registry.gaussian_not_optional(
          scatterer.scattering_type);
      FloatType fdp = scatterer.fdp;
      fractional<FloatType> coor_frac = scatterer.site;
      detail::gaussian_fourier_transformed<FloatType> gaussian_ft(
        exp_table,
        gaussian, scatterer.fp, scatterer.fdp, scatterer.weight(),
        *u_radius++, this->u_extra_);
      detail::calc_box<FloatType, grid_point_type, XrayScattererType>
        sampling_box(
          this->unit_cell_, this->rho_cutoff_, this->max_d_sq_upper_bound_,
          grid_f, coor_frac, gaussian_ft);
      if (sampling_box.excessive_radius) {
        this->excessive_sampling_radius_i_seqs_.push_back(i_seq);
      }
      this->update_sampling_box_statistics(
        sampling_box.n_points, sampling_box.box_edges);
      if (sampled_density_must_be_positive) {
        if (   gaussian_ft.rho_real_0() < 0
            || gaussian_ft.rho_real(sampling_box.max_d_sq) < 0) {

          throw error("Negative electron density at sampling point.");
        }
      }
      if (scatterer.flags.use_u_aniso()) {
        gaussian_ft = detail::gaussian_fourier_transformed<FloatType>(
          exp_table,
          gaussian, scatterer.fp, scatterer.fdp, scatterer.weight(),
          *u_cart++, this->u_extra_);
      }
      if (gifes != 0) {
        grid_indices_for_one_scatterer.reserve(sampling_box.n_points);
      }
      std::size_t exp_tab_size = exp_table.table_.size();
#     include <cctbx/xray/sampling_loop.h>
        if (gifes != 0) {
          grid_indices_for_one_scatterer.push_back(i_map);
          if (store_grid_indices_for_each_scatterer < 0) continue;
        }
        if (this->anomalous_flag_) i_map *= 2;
        if (!scatterer.flags.use_u_aniso()) {
#ifdef CCTBX_READABLE_CODE
          map_begin[i_map] += gaussian_ft.rho_real(d_sq);
#else
          if (exp_table.one_over_step_size_ == 0) {
            FloatType contr = 0;
            for (std::size_t i=0;i<gaussian_ft.n_rho_real_terms;i++) {
              contr += gaussian_ft.as_real_[i]
                     * std::exp(gaussian_ft.bs_real_[i] * d_sq);
            }
            map_begin[i_map] += contr;
          }
          else {
            FloatType contr = 0;
            FloatType d_sq_et = d_sq * exp_table.one_over_step_size_;
            for (std::size_t i=0;i<gaussian_ft.n_rho_real_terms;i++) {
              FloatType xs = gaussian_ft.bs_real_[i] * d_sq_et;
              std::size_t j = static_cast<std::size_t>(xs+.5);
              FloatType e;
              if (j >= exp_tab_size) {
                exp_table.expand(j + 1);
                exp_tab_size = exp_table.table_.size();
                e = exp_table.table_[j];
              }
              else {
                e = exp_table.table_[j];
              }
              contr += gaussian_ft.as_real_[i] * e;
            }
            /*
            Ideally an assertion i_map<real_map_.size() should be here.
            Perhaps it is not here for performance reasons.
            */
            map_begin[i_map] += contr;
          }
#endif
          if (fdp) {
            map_begin[i_map+1] += gaussian_ft.rho_imag(d_sq);
          }
        }
        else {
          map_begin[i_map] += gaussian_ft.rho_real(d);
          if (fdp) {
            map_begin[i_map+1] += gaussian_ft.rho_imag(d);
          }
        }
      CCTBX_XRAY_SAMPLING_LOOP_END
      if (gifes != 0) {
        gifes[i_seq].assign(
          grid_indices_for_one_scatterer.begin(),
          grid_indices_for_one_scatterer.end());
        grid_indices_for_one_scatterer.clear();
      }
    }
    CCTBX_ASSERT(u_radius == this->u_radius_cache_.end());
    CCTBX_ASSERT(u_cart == this->u_cart_cache_.end());
    this->u_radius_cache_ = u_radius_cache_t(); // free memory
    this->u_cart_cache_ = u_cart_cache_t(); // free memory
    this->exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H
