/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Copyright (c) 2002 Airlie McCoy.

   Revision history:
     2002 Nov: Modified fragments from cctbx/sftbx/sfmap.h (rwgk)
     2002 May: Created based on phaser/src/MapFFT.cc by Airlie McCoy (rwgk)
 */

#ifndef CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H
#define CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/maptbx/grid_tags.h>

namespace cctbx { namespace xray {

  //! Artificial temperature factor for the treatment of aliasing problems.
  /*! Reference:

        Gerard Bricogne (2001),
        International Tables for Crystallography, Volume B, 2001, p. 87
        (end of section 1.3.4.4.5).

      @param d_min = 1/d*max
      @param grid_resolution_factor = 1/(2*sigma)
      @param quality_factor = Q
      @param max_u_extra is a user-defined upper limit.

      quality_factor = 100 for 1% accuracy.

      See also:
        cctbx::uctbx::unit_cell::d()
    */
  inline // Potential NOINLINE
  double
  calc_u_extra(double d_min,
               double grid_resolution_factor,
               double quality_factor=100,
               double max_u_extra=adptbx::b_as_u(1000))
  {
    CCTBX_ASSERT(d_min > 0);
    double numerator = adptbx::b_as_u(std::log10(quality_factor));
    double sigma = 1 / (2 * grid_resolution_factor);
    double denominator = sigma * (sigma - 1) / (d_min * d_min);
    if (max_u_extra * denominator > numerator) {
      return numerator / denominator;
    }
    return max_u_extra;
  }

  template <typename FloatType>
  void
  eliminate_u_extra(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_extra,
    af::const_ref<miller::index<> > const& miller_indices,
    af::ref<std::complex<FloatType> > const& structure_factors,
    FloatType const& norm=1)
  {
    CCTBX_ASSERT(miller_indices.size() == structure_factors.size());
    //    f *= exp(b_extra * d_star_sq(h)/4)
    // => f *= exp(2*pi*pi*u_extra * d_star_sq(h))
    FloatType tppu = scitbx::constants::two_pi_sq * u_extra;
    for(std::size_t i=0;i<miller_indices.size();i++) {
      structure_factors[i] *= norm * FloatType(std::exp(
        tppu * unit_cell.d_star_sq(miller_indices[i])));
    }
  }

  namespace detail {

    // self-expanding exponent table
    template <typename FloatType>
    class exponent_table
    {
      public:
        exponent_table() {}

        explicit
        exponent_table(FloatType const& one_over_step_size)
        : one_over_step_size_(one_over_step_size)
        {}

        FloatType
        operator()(FloatType const& x)
        {
          if (one_over_step_size_ == 0) return std::exp(x);
          FloatType xs = x * one_over_step_size_;
          CCTBX_ASSERT(xs >= 0);
          std::size_t i(xs + FloatType(.5));
          if (i >= table_.size()) expand(i + 1);
          return table_[i];
        }

        std::vector<FloatType> const&
        table() const { return table_; }

      private:
        FloatType one_over_step_size_;
        std::vector<FloatType> table_;

        void expand(std::size_t n);
    };

    template <typename FloatType>
    void
    exponent_table<FloatType>::expand(std::size_t n)
    {
      table_.reserve(n);
      for(std::size_t i = table_.size(); i < n; i++) {
        table_.push_back(std::exp(i / one_over_step_size_));
      }
    }

    // ff(hc) =
    //   a Exp[-hc.b_arg.hc]
    // rho_anisotropic(rc) =
    //   a / (2 Sqrt[2] Sqrt[Det[b_arg]]) Exp[-1/4 rc.Inverse[b_arg].rc]
    // XXX explanation for additional factors?
    template <typename FloatType>
    inline
    void
    anisotropic_3d_gaussian_fourier_transform(
      FloatType const& a,
      scitbx::sym_mat3<FloatType> const& b_all,
      FloatType& as,
      scitbx::sym_mat3<FloatType>& bs)
    {
      using scitbx::constants::pi;
      using scitbx::constants::pi_sq;
      FloatType d = b_all.determinant();
      CCTBX_ASSERT(d != 0);
      scitbx::sym_mat3<FloatType> cfmt = b_all.co_factor_matrix_transposed();
      as = a / ((2 * std::sqrt(2.)) * std::sqrt(d))
         * (16 * std::sqrt(2.) * pi * std::sqrt(pi));
      bs = cfmt / FloatType(-4 * d) * FloatType(16 * pi_sq);
    }

    template <typename FloatType>
    inline
    void
    isotropic_3d_gaussian_fourier_transform(
      FloatType const& a,
      FloatType const& b_all,
      FloatType& as,
      FloatType& bs)
    {
      using scitbx::constants::pi;
      using scitbx::constants::pi_sq;
      FloatType d = b_all * b_all * b_all;
      CCTBX_ASSERT(d != 0);
      as = a / ((2 * std::sqrt(2.)) * std::sqrt(d))
         * (16 * std::sqrt(2.) * pi * std::sqrt(pi));
      bs = 1 / (-4 * b_all) * FloatType(16 * pi_sq);
    }

    template <typename FloatTypeCaasfB,
              typename FloatType>
    inline
    scitbx::sym_mat3<FloatType>
    compose_anisotropic_b_all(
      FloatTypeCaasfB const& caasf_b,
      FloatType const& u_extra,
      scitbx::sym_mat3<FloatType> const& u_cart)
    {
      return scitbx::sym_mat3<FloatType>(caasf_b + adptbx::u_as_b(u_extra))
           + scitbx::sym_mat3<FloatType>(adptbx::u_as_b(u_cart));
    }

    template <typename FloatType,
              typename CaasfType>
    class caasf_fourier_transformed
    {
      public:
        caasf_fourier_transformed() {}

        caasf_fourier_transformed(
          CaasfType const& caasf,
          std::complex<FloatType> const& fp_fdp,
          FloatType const& w,
          FloatType const& u_iso,
          FloatType const& u_extra)
          : anisotropic_flag_(false)
        {
          FloatType b_incl_extra = adptbx::u_as_b(u_iso + u_extra);
          std::size_t i = 0;
          for(;i<caasf.n_ab();i++) {
            isotropic_3d_gaussian_fourier_transform(
              w * caasf.a(i), caasf.b(i) + b_incl_extra,
              as_real_[i], bs_real_[i]);
          }
          isotropic_3d_gaussian_fourier_transform(
            w * (caasf.c() + fp_fdp.real()), b_incl_extra,
            as_real_[i], bs_real_[i]);
          isotropic_3d_gaussian_fourier_transform(
            w * (fp_fdp.imag()), b_incl_extra,
            as_imag_, bs_imag_);
        }

        caasf_fourier_transformed(
          CaasfType const& caasf,
          std::complex<FloatType> const& fp_fdp,
          FloatType const& w,
          scitbx::sym_mat3<FloatType> const& u_cart,
          FloatType const& u_extra)
          : anisotropic_flag_(true)
        {
          std::size_t i = 0;
          for(;i<caasf.n_ab();i++) {
            anisotropic_3d_gaussian_fourier_transform(
              w * caasf.a(i),
              compose_anisotropic_b_all(caasf.b(i), u_extra, u_cart),
              as_real_[i], aniso_bs_real_[i]);
          }
          anisotropic_3d_gaussian_fourier_transform(
            w * (caasf.c() + fp_fdp.real()),
            compose_anisotropic_b_all(0, u_extra, u_cart),
            as_real_[i], aniso_bs_real_[i]);
          anisotropic_3d_gaussian_fourier_transform(
            w * (fp_fdp.imag()),
            compose_anisotropic_b_all(0, u_extra, u_cart),
            as_imag_, aniso_bs_imag_);
        }

        bool
        anisotropic_flag() const { return anisotropic_flag_; }

        FloatType
        rho_real_0() const
        {
          return af::sum(as_real_);
        }

        FloatType
        rho_real(exponent_table<FloatType>& exp_table,
                 FloatType const& d_sq) const
        {
          FloatType r(0);
          for (std::size_t i=0;i<as_real_.size();i++) {
            r += as_real_[i] * exp_table(bs_real_[i] * d_sq);
          }
          return r;
        }

        FloatType
        rho_real(exponent_table<FloatType>& exp_table,
                 scitbx::vec3<FloatType> const& d) const
        {
          FloatType r(0);
          for (std::size_t i=0;i<as_real_.size();i++) {
            r += as_real_[i] * exp_table(d * aniso_bs_real_[i] * d);
          }
          return r;
        }

        FloatType
        rho_imag(exponent_table<FloatType>& exp_table,
                 FloatType const& d_sq) const
        {
          return as_imag_ * exp_table(bs_imag_ * d_sq);
        }

        FloatType
        rho_imag(exponent_table<FloatType>& exp_table,
                 scitbx::vec3<FloatType> const& d) const
        {
          return as_imag_ * exp_table(d * aniso_bs_imag_ * d);
        }

      protected:
        bool anisotropic_flag_;
        af::tiny<FloatType, CaasfType::n_plus_1> as_real_;
        af::tiny<FloatType, CaasfType::n_plus_1> bs_real_;
        af::tiny<scitbx::sym_mat3<FloatType>, CaasfType::n_plus_1>
          aniso_bs_real_;
        FloatType as_imag_;
        FloatType bs_imag_;
        scitbx::sym_mat3<FloatType> aniso_bs_imag_;
    };

    template <typename FloatType,
              typename GridPointType>
    struct calc_shell
    {
      GridPointType radii;
      FloatType max_d_sq;

      template <typename CaasfType>
      calc_shell(
        uctbx::unit_cell const& unit_cell,
        FloatType const& wing_cutoff,
        GridPointType const& grid_n,
        caasf_fourier_transformed<FloatType, CaasfType> const& caasf_ft,
        exponent_table<FloatType>& exp_table)
      :
        max_d_sq(0)
      {
        CCTBX_ASSERT(!caasf_ft.anisotropic_flag());
        af::tiny<FloatType, 3> grid_n_f = grid_n;
        FloatType rho_cutoff = scitbx::fn::absolute(
          caasf_ft.rho_real_0() * wing_cutoff);
        for(std::size_t i_basis_vec=0;i_basis_vec<3;i_basis_vec++) {
          typename GridPointType::value_type ig;
          for (ig = 1; ig < grid_n[i_basis_vec]; ig++) {
            fractional<FloatType> d_frac(0,0,0);
            d_frac[i_basis_vec] = ig / grid_n_f[i_basis_vec];
            FloatType d_sq = unit_cell.length_sq(d_frac);
            if (scitbx::fn::absolute(
                  caasf_ft.rho_real(exp_table, d_sq)) < rho_cutoff) {
              break;
            }
            if (max_d_sq < d_sq) max_d_sq = d_sq;
          }
          if (ig == 1) {
            throw error(
              "Error determining sampling radius for"
              " model electron density. This may be due to incorrectly"
              " defined atomic scattering factors.");
          }
          if (ig == grid_n[i_basis_vec]) {
            throw error(
              "Excessive radius for real-space sampling of"
              " model electron density.");
          }
          radii[i_basis_vec] = ig - 1;
        }
      }
    };

    template <typename FloatType>
    inline
    FloatType
    add_for_rounding(FloatType const& x)
    {
      if (x < FloatType(0)) return x - FloatType(.5);
      return x + FloatType(.5);
    }

    template <typename FloatType,
              typename GridPointType>
    inline
    GridPointType
    calc_nearest_grid_point(fractional<FloatType> const& coor,
                            GridPointType const& grid_n)
    {
      typedef typename GridPointType::value_type value_type;
      GridPointType grid_point;
      for(std::size_t i=0;i<3;i++) {
        grid_point[i] = value_type(add_for_rounding(coor[i] * grid_n[i]));
      }
      return grid_point;
    }

    template <typename ArrayMaxType,
              typename ArrayValType>
    ArrayMaxType&
    array_update_max(ArrayMaxType& am, ArrayValType const& av)
    {
      for(std::size_t i=0;i<am.size();i++) {
        if (am[i] < av[i]) am[i] = av[i];
      }
      return am;
    }

  } // namespace detail

  template <typename FloatType = double,
            typename XrayScattererType = scatterer<> >
  class sampled_model_density
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
    grid_point_type const& grid_f = map_accessor_.focus();
    grid_point_type const& grid_a = map_accessor_.all();
    detail::exponent_table<FloatType> exp_table(
      exp_table_one_over_step_size);
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
      }}}
    }
    exp_table_size_ = exp_table.table().size();
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SAMPLED_MODEL_DENSITY_H
