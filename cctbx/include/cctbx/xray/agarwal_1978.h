#ifndef CCTBX_XRAY_AGARWAL_1978_H
#define CCTBX_XRAY_AGARWAL_1978_H

#include <cctbx/xray/sampled_model_density.h>

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
          caasf_fourier_transformed<FloatType, CaasfType> const& base)
        :
          caasf_fourier_transformed<FloatType, CaasfType>(base)
        {}

        scitbx::vec3<FloatType>
        d_rho_real_d_site(exponent_table<FloatType>& exp_table,
                          scitbx::vec3<FloatType> const& d,
                          FloatType const& d_sq) const
        {
          scitbx::vec3<FloatType> drds(0,0,0);
          for (std::size_t i=0;i<this->as_real_.size();i++) {
            drds += d * (2 * this->bs_real_[i]
                         * this->as_real_[i]
                         * exp_table(this->bs_real_[i] * d_sq));
          }
          return drds;
        }
    };

  } // namespace detail

  template <typename FloatType = double,
            typename XrayScattererType = scatterer<> >
  class agarwal_1978
  {
    public:
      typedef XrayScattererType xray_scatterer_type;
      typedef typename xray_scatterer_type::caasf_type::base_type caasf_type;

      typedef typename maptbx::c_grid_padded_p1<3> accessor_type;
      typedef typename accessor_type::index_type grid_point_type;
      typedef typename grid_point_type::value_type grid_point_element_type;

      typedef af::versa<FloatType, accessor_type>
        real_map_type;
      typedef af::versa<std::complex<FloatType>, accessor_type>
        complex_map_type;

      agarwal_1978() {}

      agarwal_1978(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        af::const_ref<std::complex<FloatType>, accessor_type> const& ft_dt,
        bool lifchitz,
        grid_point_type const& fft_n_real,
        grid_point_type const& fft_m_real,
        FloatType const& u_extra=0.25,
        FloatType const& wing_cutoff=1.e-3,
        FloatType const& exp_table_one_over_step_size=-100,
        bool force_complex=false,
        bool electron_density_must_be_positive=true);

      uctbx::unit_cell const&
      unit_cell() { return unit_cell_; }

      FloatType
      u_extra() const { return u_extra_; }

      FloatType
      wing_cutoff() const { return wing_cutoff_; }

      FloatType
      exp_table_one_over_step_size() const
      {
        return exp_table_one_over_step_size_;
      }

      std::size_t
      n_scatterers_passed() const { return n_scatterers_passed_; }

      std::size_t
      n_contributing_scatterers() const { return n_contributing_scatterers_; }

      std::size_t
      n_anomalous_scatterers() const { return n_anomalous_scatterers_; }

      bool
      anomalous_flag() const { return complex_map_.size() != 0; }

      real_map_type
      real_map() { return real_map_; }

      complex_map_type
      complex_map() { return complex_map_; }

      std::size_t
      exp_table_size() const { return exp_table_size_; }

      grid_point_type const&
      max_shell_radii() const { return max_shell_radii_; }

      fractional<FloatType>
      max_shell_radii_frac() const
      {
        fractional<FloatType> r;
        for(std::size_t i=0;i<3;i++) {
          r[i] = FloatType(max_shell_radii_[i]) / map_accessor_.focus()[i];
        }
        return r;
      }

      af::shared<std::complex<double> >
      grad() const { return grad_; }

      af::shared<std::complex<double> >
      grad_x() const { return grad_x_; }

      af::shared<std::complex<double> >
      grad_y() const { return grad_y_; }

      af::shared<std::complex<double> >
      grad_z() const { return grad_z_; }

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
        FloatType norm = unit_cell_.volume() / map_accessor_.focus_size_1d();
        eliminate_u_extra(
          unit_cell_, u_extra_, miller_indices, structure_factors, norm);
      }

    private:
      uctbx::unit_cell unit_cell_;
      std::size_t n_scatterers_passed_;
      FloatType u_extra_;
      FloatType wing_cutoff_;
      FloatType exp_table_one_over_step_size_;
      std::size_t n_anomalous_scatterers_;
      std::size_t n_contributing_scatterers_;
      accessor_type map_accessor_;
      real_map_type real_map_;
      complex_map_type complex_map_;
      std::size_t exp_table_size_;
      grid_point_type max_shell_radii_;
      af::shared<std::complex<double> > grad_;
      af::shared<std::complex<double> > grad_x_;
      af::shared<std::complex<double> > grad_y_;
      af::shared<std::complex<double> > grad_z_;
  };

  template <typename FloatType,
            typename XrayScattererType>
  agarwal_1978<FloatType, XrayScattererType>
  ::agarwal_1978(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<XrayScattererType> const& scatterers,
    af::const_ref<std::complex<FloatType>, accessor_type> const& ft_dt,
    bool lifchitz,
    grid_point_type const& fft_n_real,
    grid_point_type const& fft_m_real,
    FloatType const& u_extra,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    bool force_complex,
    bool electron_density_must_be_positive)
  :
    unit_cell_(unit_cell),
    n_scatterers_passed_(scatterers.size()),
    u_extra_(u_extra),
    wing_cutoff_(wing_cutoff),
    exp_table_one_over_step_size_(exp_table_one_over_step_size),
    n_contributing_scatterers_(0),
    n_anomalous_scatterers_(0),
    exp_table_size_(0),
    max_shell_radii_(0,0,0)
  {
    scitbx::mat3<FloatType> orth_mx = unit_cell_.orthogonalization_matrix();
    if (orth_mx[3] != 0 || orth_mx[6] != 0 || orth_mx[7] != 0) {
      throw error(
        "Fatal Programming Error:"
        " Real-space sampling of model electron density"
        " is optimized for orthogonalization matrix"
        " according to the PDB convention. The orthogonalization"
        " matrix passed is not compatible with this convention.");
    }
    const xray_scatterer_type* scatterer;
    for(scatterer=scatterers.begin();scatterer!=scatterers.end();scatterer++) {
      if (scatterer->weight() == 0) continue;
      n_contributing_scatterers_++;
      if (scatterer->fp_fdp.imag() != 0) {
        n_anomalous_scatterers_++;
      }
    }
    bool anomalous_flag;
    FloatType* map_begin;
    if (n_anomalous_scatterers_ == 0 && !force_complex) {
      anomalous_flag = false;
      map_accessor_ = accessor_type(fft_m_real, fft_n_real);
      real_map_.resize(map_accessor_);
      map_begin = real_map_.begin();
    }
    else {
      anomalous_flag = true;
      map_accessor_ = accessor_type(fft_n_real, fft_n_real);
      complex_map_.resize(map_accessor_);
      map_begin = reinterpret_cast<FloatType*>(complex_map_.begin());
    }
    CCTBX_ASSERT(anomalous_flag == false); // XXX
    grid_point_type const& grid_f = map_accessor_.focus();
    grid_point_type const& grid_a = map_accessor_.all();
    detail::exponent_table<FloatType> exp_table(
      exp_table_one_over_step_size);
    for(scatterer=scatterers.begin();scatterer!=scatterers.end();scatterer++)
    {
      CCTBX_ASSERT(!scatterer->anisotropic_flag); // XXX
      if (scatterer->weight() == 0) continue;
      FloatType fdp = scatterer->fp_fdp.imag();
      fractional<FloatType> coor_frac = scatterer->site;
      FloatType u_iso;
      scitbx::sym_mat3<FloatType> u_cart;
      if (!scatterer->anisotropic_flag) {
        u_iso = scatterer->u_iso;
      }
      else {
        u_cart = adptbx::u_star_as_u_cart(unit_cell_, scatterer->u_star);
        scitbx::vec3<FloatType> ev = adptbx::eigenvalues(u_cart);
        CCTBX_ASSERT(adptbx::is_positive_definite(ev));
        u_iso = af::max(ev);
      }
      CCTBX_ASSERT(u_iso >= 0);
      detail::caasf_fourier_transformed<FloatType, caasf_type> caasf_ft(
        scatterer->caasf, scatterer->fp_fdp, scatterer->weight(),
        u_iso, u_extra_);
      detail::calc_shell<FloatType, grid_point_type> shell(
        unit_cell_, wing_cutoff_, grid_f, caasf_ft, exp_table);
      detail::array_update_max(max_shell_radii_, shell.radii);
      if (electron_density_must_be_positive) {
        if (   caasf_ft.rho_real_0() < 0
            || caasf_ft.rho_real(exp_table, shell.max_d_sq) < 0) {

          throw error("Negative electron density at sampling point.");
        }
      }
      if (scatterer->anisotropic_flag) {
        caasf_ft = detail::caasf_fourier_transformed<FloatType,caasf_type>(
          scatterer->caasf, scatterer->fp_fdp, scatterer->weight(),
          u_cart, u_extra_);
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
      std::complex<FloatType> gr(0);
      scitbx::vec3<std::complex<FloatType> > grs(0,0,0);
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
        if (anomalous_flag) i_map *= 2;
        if (!scatterer->anisotropic_flag) {
          map_begin[i_map] += caasf_ft.rho_real(exp_table, d_sq);
          if (fdp) {
            map_begin[i_map+1] += caasf_ft.rho_imag(exp_table, d_sq);
          }
        }
        else {
          map_begin[i_map] += caasf_ft.rho_real(exp_table, d);
          if (fdp) {
            map_begin[i_map+1] += caasf_ft.rho_imag(exp_table, d);
          }
        }
        if (!lifchitz) {
          gr += ft_dt(gp) * map_begin[i_map];
        }
        else {
          scitbx::vec3<FloatType>
            drds = detail::d_caasf_fourier_transformed<FloatType,caasf_type>(
              caasf_ft).d_rho_real_d_site(exp_table, d, d_sq)
            * unit_cell_.orthogonalization_matrix();
          for(std::size_t i=0;i<3;i++) {
            grs[i] += ft_dt(gp) * drds[i];
          }
        }
      }}}
      if (!lifchitz) {
        grad_.push_back(gr);
      }
      else {
        grad_x_.push_back(grs[0]);
        grad_y_.push_back(grs[1]);
        grad_z_.push_back(grs[2]);
      }
    }
    exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_AGARWAL_1978_H
