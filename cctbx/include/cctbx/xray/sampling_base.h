#ifndef CCTBX_XRAY_SAMPLING_BASE_H
#define CCTBX_XRAY_SAMPLING_BASE_H

#include <cctbx/xray/scatterer.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>

namespace cctbx { namespace xray {

  namespace detail {

    using scitbx::constants::four_pi_sq;
    using scitbx::constants::eight_pi_sq;
    static const double eight_pi_pow_3_2
      = 8 * scitbx::constants::pi * std::sqrt(scitbx::constants::pi);

  } // namespace detail

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

  // Not available in Python.
  template <typename FloatType>
  inline
  void
  apply_u_extra(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_extra,
    miller::index<> const& miller_index,
    std::complex<FloatType>& structure_factor,
    FloatType const& multiplier)
  {
    //    f *= exp(b_extra * d_star_sq(h)/4)
    // => f *= exp(2*pi*pi*u_extra * d_star_sq(h))
    FloatType tppu = scitbx::constants::two_pi_sq * u_extra;
    structure_factor *= multiplier * FloatType(std::exp(
        tppu * unit_cell.d_star_sq(miller_index)));
  }

  template <typename FloatType>
  void
  apply_u_extra(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_extra,
    af::const_ref<miller::index<> > const& miller_indices,
    af::ref<std::complex<FloatType> > const& structure_factors,
    FloatType const& multiplier=1)
  {
    CCTBX_ASSERT(miller_indices.size() == structure_factors.size());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      apply_u_extra(
        unit_cell, u_extra, miller_indices[i], structure_factors[i],
        multiplier);
    }
  }

  template <typename FloatType>
  void
  apply_u_extra(
    uctbx::unit_cell const& unit_cell,
    FloatType const& u_extra,
    af::const_ref<miller::index<> > const& miller_indices,
    af::ref<std::complex<FloatType> > const& structure_factors,
    af::const_ref<FloatType> const& multipliers)
  {
    CCTBX_ASSERT(miller_indices.size() == structure_factors.size());
    CCTBX_ASSERT(miller_indices.size() == multipliers.size());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      apply_u_extra(
        unit_cell, u_extra, miller_indices[i], structure_factors[i],
        multipliers[i]);
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

    // ff(hc) = a Exp[-hc.b_arg.hc]
    // rho_anisotropic(rc) =
    //   8 Pi^(3/2) a / Sqrt[Det[b_arg]] Exp[-4 Pi^2 rc.Inverse[b_arg].rc]
    template <typename FloatType>
    inline
    void
    anisotropic_3d_gaussian_fourier_transform(
      FloatType const& a,
      scitbx::sym_mat3<FloatType> const& b_all,
      FloatType& as,
      scitbx::sym_mat3<FloatType>& bs)
    {
      FloatType d = b_all.determinant();
      CCTBX_ASSERT(d != 0);
      scitbx::sym_mat3<FloatType> cfmt = b_all.co_factor_matrix_transposed();
      as = eight_pi_pow_3_2 * a / std::sqrt(d);
      bs = cfmt / (d / -four_pi_sq);
    }

    template <typename FloatType>
    inline
    FloatType
    anisotropic_3d_gaussian_fourier_transform(
      FloatType const& a,
      scitbx::sym_mat3<FloatType> const& b_all)
    {
      FloatType d = b_all.determinant();
      CCTBX_ASSERT(d != 0);
      return eight_pi_pow_3_2 * a / std::sqrt(d);
    }

    template <typename FloatType>
    inline
    FloatType
    isotropic_3d_gaussian_fourier_transform(
      FloatType const& a,
      FloatType const& b_all)
    {
      FloatType d = b_all * b_all * b_all;
      CCTBX_ASSERT(d != 0);
      return eight_pi_pow_3_2 * a / std::sqrt(d);
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
      as = isotropic_3d_gaussian_fourier_transform(a, b_all);
      bs = -four_pi_sq / b_all;
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
        BOOST_STATIC_CONSTANT(
          std::size_t, n_rho_real_terms = CaasfType::n_plus_1);

        caasf_fourier_transformed() {}

        caasf_fourier_transformed(
          exponent_table<FloatType>& exp_table,
          CaasfType const& caasf,
          std::complex<FloatType> const& fp_fdp,
          FloatType const& w,
          FloatType const& u_iso,
          FloatType const& u_extra)
        :
          exp_table_(&exp_table),
          anisotropic_flag_(false)
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
          exponent_table<FloatType>& exp_table,
          CaasfType const& caasf,
          std::complex<FloatType> const& fp_fdp,
          FloatType const& w,
          scitbx::sym_mat3<FloatType> const& u_cart,
          FloatType const& u_extra)
        :
          exp_table_(&exp_table),
          anisotropic_flag_(true)
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

        FloatType
        exp_table(FloatType const& x) const
        {
          return (*exp_table_)(x);
        }

        bool
        anisotropic_flag() const { return anisotropic_flag_; }

        FloatType
        rho_real_0() const
        {
          return af::sum(as_real_);
        }

        FloatType
        exp_term(FloatType const& d_sq, std::size_t i) const
        {
          return exp_table(bs_real_[i] * d_sq);
        }

        FloatType
        exp_term(scitbx::vec3<FloatType> const& d, std::size_t i) const
        {
          return exp_table(d * aniso_bs_real_[i] * d);
        }

        FloatType
        exp_term(FloatType const& d_sq) const
        {
          return exp_table(bs_imag_ * d_sq);
        }

        FloatType
        exp_term(scitbx::vec3<FloatType> const& d) const
        {
          return exp_table(d * aniso_bs_imag_ * d);
        }

        template <typename DistanceType>
        FloatType
        rho_real_term(DistanceType const& d_or_d_sq, std::size_t i) const
        {
          return as_real_[i] * exp_term(d_or_d_sq, i);
        }

        template <typename DistanceType>
        FloatType
        rho_real(DistanceType const& d_or_d_sq) const
        {
          FloatType r(0);
          for (std::size_t i=0;i<n_rho_real_terms;i++) {
            r += rho_real_term(d_or_d_sq, i);
          }
          return r;
        }

        template <typename DistanceType>
        FloatType
        rho_imag(DistanceType const& d_or_d_sq) const
        {
          return as_imag_ * exp_term(d_or_d_sq);
        }

      protected:
        exponent_table<FloatType>* exp_table_;
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
        caasf_fourier_transformed<FloatType, CaasfType> const& caasf_ft)
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
                  caasf_ft.rho_real(d_sq)) < rho_cutoff) {
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

  template <typename FloatType=double,
            typename XrayScattererType=scatterer<> >
  class sampling_base
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

      sampling_base() {}

      sampling_base(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        FloatType const& u_extra,
        FloatType const& wing_cutoff,
        FloatType const& exp_table_one_over_step_size,
        FloatType const& tolerance_positive_definite);

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

      FloatType
      tolerance_positive_definite() const
      {
        return tolerance_positive_definite_;
      }

      std::size_t
      n_scatterers_passed() const { return n_scatterers_passed_; }

      std::size_t
      n_contributing_scatterers() const { return n_contributing_scatterers_; }

      std::size_t
      n_anomalous_scatterers() const { return n_anomalous_scatterers_; }

      bool
      anomalous_flag() const { return anomalous_flag_; }

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

    protected:
      uctbx::unit_cell unit_cell_;
      std::size_t n_scatterers_passed_;
      FloatType u_extra_;
      FloatType wing_cutoff_;
      FloatType exp_table_one_over_step_size_;
      FloatType tolerance_positive_definite_;
      std::size_t n_contributing_scatterers_;
      std::size_t n_anomalous_scatterers_;
      bool anomalous_flag_;
      std::size_t exp_table_size_;
      grid_point_type max_shell_radii_;
      accessor_type map_accessor_;

      FloatType
      get_u_cart_and_u_iso(
        scitbx::sym_mat3<FloatType> const& u_star,
        scitbx::sym_mat3<FloatType>& u_cart)
      {
        u_cart = adptbx::u_star_as_u_cart(unit_cell_, u_star);
        scitbx::vec3<FloatType>
          u_cart_eigenvalues = adptbx::eigenvalues(u_cart);
        CCTBX_ASSERT(adptbx::is_positive_definite(u_cart_eigenvalues,
          tolerance_positive_definite_));
        return af::max(u_cart_eigenvalues);
      }
  };

  template <typename FloatType,
            typename XrayScattererType>
  sampling_base<FloatType, XrayScattererType>
  ::sampling_base(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<XrayScattererType> const& scatterers,
    FloatType const& u_extra,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    FloatType const& tolerance_positive_definite)
  :
    unit_cell_(unit_cell),
    n_scatterers_passed_(scatterers.size()),
    u_extra_(u_extra),
    wing_cutoff_(wing_cutoff),
    exp_table_one_over_step_size_(exp_table_one_over_step_size),
    tolerance_positive_definite_(tolerance_positive_definite),
    n_contributing_scatterers_(0),
    n_anomalous_scatterers_(0),
    anomalous_flag_(false),
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
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SAMPLING_BASE_H
