#ifndef CCTBX_XRAY_SAMPLING_BASE_H
#define CCTBX_XRAY_SAMPLING_BASE_H

#include <cctbx/xray/scattering_type_registry.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>

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
      @param max_u_base is a user-defined upper limit.

      quality_factor = 100 for 1% accuracy.

      See also:
        cctbx::uctbx::unit_cell::d()
    */
  inline // Potential NOINLINE
  double
  calc_u_base(double d_min,
              double grid_resolution_factor,
              double quality_factor=100,
              double max_u_base=adptbx::b_as_u(1000))
  {
    CCTBX_ASSERT(d_min > 0);
    double numerator = adptbx::b_as_u(std::log10(quality_factor));
    double sigma = 1 / (2 * grid_resolution_factor);
    double denominator = sigma * (sigma - 1) / (d_min * d_min);
    if (max_u_base * denominator > numerator) {
      return numerator / denominator;
    }
    return max_u_base;
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

  namespace detail {

    // self-expanding exponent table
    template <typename FloatType>
    class exponent_table
    {
      public:
        static const std::size_t excessive_range_error_limit = 1000000;

        exponent_table() {}

        explicit
        exponent_table(
          FloatType const& one_over_step_size,
          std::size_t initial_reserve=10000)//avoids reallocation in most cases
        :
          one_over_step_size_(one_over_step_size)
        {
          if (one_over_step_size_ != 0) {
            table_.reserve(initial_reserve);
          }
        }

        FloatType
        operator()(FloatType const& x)
        {
          if (one_over_step_size_ == 0) return std::exp(x);
          FloatType xs = x * one_over_step_size_;
          std::size_t i = static_cast<std::size_t>(xs+.5);
          if (i >= table_.size()) {
            expand(i + 1);
          }
          return table_[i];
        }

        std::vector<FloatType> const&
        table() const { return table_; }

      public: // not protected to allow for manually inlined code
        FloatType one_over_step_size_;
        std::vector<FloatType> table_;

        void expand(std::size_t n);
    };

    template <typename FloatType>
    void
    exponent_table<FloatType>::expand(std::size_t n)
    {
      if (n > excessive_range_error_limit) {
        throw std::runtime_error(
          __FILE__ ": exponent_table: excessive range.");
      }
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
      CCTBX_ASSERT(d > 0);
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
      CCTBX_ASSERT(d > 0);
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
      if (d <= 0) {
        char buf[80];
        std::snprintf(buf, sizeof(buf),
          "isotropic_3d_gaussian_fourier_transform: b_all=%.6g", b_all);
        throw CCTBX_ERROR(buf);
      }
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

    template <typename FloatTypeGaussianB,
              typename FloatType>
    inline
    scitbx::sym_mat3<FloatType>
    compose_anisotropic_b_all(
      FloatTypeGaussianB const& gaussian_b,
      FloatType const& u_extra,
      scitbx::sym_mat3<FloatType> const& u_cart)
    {
      return scitbx::sym_mat3<FloatType>(gaussian_b + adptbx::u_as_b(u_extra))
           + scitbx::sym_mat3<FloatType>(adptbx::u_as_b(u_cart));
    }

    template <typename FloatType>
    class gaussian_fourier_transformed
    {
      public:
        BOOST_STATIC_CONSTANT(std::size_t, max_n_rho_real_terms =
          eltbx::xray_scattering::gaussian::max_n_terms+1);

        gaussian_fourier_transformed() {}

        gaussian_fourier_transformed(
          exponent_table<FloatType>& exp_table,
          eltbx::xray_scattering::gaussian const& gaussian,
          FloatType const& fp,
          FloatType const& fdp,
          FloatType const& w,
          FloatType const& u_iso,
          FloatType const& u_extra,
          bool gradient_calculations=false)
        :
          exp_table_(&exp_table),
          anisotropic_flag_(false),
          n_rho_real_terms(gaussian.n_terms())
        {
          FloatType b_incl_extra = adptbx::u_as_b(u_iso + u_extra);
          std::size_t i = 0;
          for(;i<gaussian.n_terms();i++) {
            scitbx::math::gaussian::term<double> ti = gaussian.terms()[i];
            isotropic_3d_gaussian_fourier_transform(
              w * ti.a, ti.b + b_incl_extra,
              as_real_[i], bs_real_[i]);
          }
          FloatType c_fp = gaussian.c() + fp;
          if (c_fp != 0 || gradient_calculations) {
            isotropic_3d_gaussian_fourier_transform(
              w * c_fp, b_incl_extra,
              as_real_[i], bs_real_[i]);
            n_rho_real_terms++;
          }
          if (fdp != 0 || gradient_calculations) {
            isotropic_3d_gaussian_fourier_transform(
              w * fdp, b_incl_extra,
              as_imag_, bs_imag_);
          }
          else {
            as_imag_ = 0;
          }
        }

        gaussian_fourier_transformed(
          exponent_table<FloatType>& exp_table,
          eltbx::xray_scattering::gaussian const& gaussian,
          FloatType const& fp,
          FloatType const& fdp,
          FloatType const& w,
          scitbx::sym_mat3<FloatType> const& u_cart,
          FloatType const& u_extra,
          bool gradient_calculations=false)
        :
          exp_table_(&exp_table),
          anisotropic_flag_(true),
          n_rho_real_terms(gaussian.n_terms())
        {
          std::size_t i = 0;
          for(;i<gaussian.n_terms();i++) {
            scitbx::math::gaussian::term<double> ti = gaussian.terms()[i];
            anisotropic_3d_gaussian_fourier_transform(
              w * ti.a, compose_anisotropic_b_all(ti.b, u_extra, u_cart),
              as_real_[i], aniso_bs_real_[i]);
          }
          FloatType c_fp = gaussian.c() + fp;
          if (c_fp != 0 || gradient_calculations) {
            anisotropic_3d_gaussian_fourier_transform(
              w * c_fp,
              compose_anisotropic_b_all(0, u_extra, u_cart),
              as_real_[i], aniso_bs_real_[i]);
            n_rho_real_terms++;
          }
          if (fdp != 0 || gradient_calculations) {
            anisotropic_3d_gaussian_fourier_transform(
              w * fdp,
              compose_anisotropic_b_all(0, u_extra, u_cart),
              as_imag_, aniso_bs_imag_);
          }
          else {
            as_imag_ = 0;
          }
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
          FloatType result = 0;
          for (std::size_t i=0;i<n_rho_real_terms;i++) {
            result += as_real_[i];
          }
          return result;
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
          if (as_imag_ == 0) return 0;
          return as_imag_ * exp_term(d_or_d_sq);
        }

        FloatType
        max_d_sq_estimate(
          FloatType const& rho_cutoff,
          FloatType const& epsilon=1.e-3) const
        {
          /* Solution of rho_cutoff = a * exp(b*max_d_sq)
             => max_d_sq = log(rho_cutoff/a) / b
           */
          if (n_rho_real_terms == 0) return 0;
          FloatType max_a = 0;
          for (std::size_t i=0;i<n_rho_real_terms;i++) {
            scitbx::math::update_max(max_a, scitbx::fn::absolute(as_real_[i]));
          }
          if (max_a <= rho_cutoff * epsilon) return 0;
          FloatType result = 0;
          for (std::size_t i=0;i<n_rho_real_terms;i++) {
            FloatType a = scitbx::fn::absolute(as_real_[i]);
            if (a <= rho_cutoff * epsilon) continue;
            scitbx::math::update_max(result,
              std::log(rho_cutoff/a)/bs_real_[i]);
          }
          return result;
        }

      public: // not protected to allow for manually inlined code
        exponent_table<FloatType>* exp_table_;
        bool anisotropic_flag_;
        std::size_t n_rho_real_terms;
        af::tiny<FloatType, max_n_rho_real_terms> as_real_;
        af::tiny<FloatType, max_n_rho_real_terms> bs_real_;
        af::tiny<scitbx::sym_mat3<FloatType>, max_n_rho_real_terms>
          aniso_bs_real_;
        FloatType as_imag_;
        FloatType bs_imag_;
        scitbx::sym_mat3<FloatType> aniso_bs_imag_;
    };

    template <typename FloatType,
              typename GridPointType,
              typename XrayScattererType>
    struct calc_box
    {
      FloatType max_d_sq;
      std::size_t n_points;
      GridPointType box_min;
      GridPointType box_max;
      GridPointType box_edges;
      bool excessive_radius;

      calc_box(
        uctbx::unit_cell const& unit_cell,
        FloatType const& rho_cutoff,
        FloatType const& max_d_sq_upper_bound,
        GridPointType const& grid_n,
        fractional<FloatType> const& coor_frac,
        gaussian_fourier_transformed<FloatType> const& gaussian_ft)
      :
        max_d_sq(0),
        n_points(1),
        excessive_radius(false)
      {
        CCTBX_ASSERT(!gaussian_ft.anisotropic_flag());
        af::tiny<FloatType, 3> grid_n_f = grid_n;
        FloatType max_d_sq_estimate = gaussian_ft.max_d_sq_estimate(
          rho_cutoff);
        for(std::size_t i_bv=0;i_bv<3;i_bv++) {
          fractional<FloatType> d_frac(0,0,0);
          d_frac[i_bv] = 1 / grid_n_f[i_bv];
          FloatType d_sq = unit_cell.length_sq(d_frac);
          FloatType gu_radius = std::sqrt(max_d_sq_estimate / d_sq);
          FloatType gu_site = coor_frac[i_bv] * grid_n_f[i_bv];
          // ensure sampling box size is at least 2x2x2
          int box_min_bound = scitbx::math::iround(gu_site);
          int box_max_bound = box_min_bound;
          if (box_max_bound < gu_site) box_max_bound++;
          else                         box_min_bound--;
          box_min[i_bv] = adjust_box_limit(
            unit_cell, rho_cutoff, max_d_sq_upper_bound,
            grid_n_f, coor_frac, gaussian_ft,
            i_bv, -1, box_min_bound,
            std::min(box_min_bound, scitbx::math::ifloor(gu_site-gu_radius)));
          box_max[i_bv] = adjust_box_limit(
            unit_cell, rho_cutoff, max_d_sq_upper_bound,
            grid_n_f, coor_frac, gaussian_ft,
            i_bv, 1, box_max_bound,
            std::max(box_max_bound, scitbx::math::iceil(gu_site+gu_radius)));
          box_edges[i_bv] = box_max[i_bv] - box_min[i_bv] + 1;
          n_points *= box_edges[i_bv];
        }
      }

      int
      adjust_box_limit(
        uctbx::unit_cell const& unit_cell,
        FloatType const& rho_cutoff,
        FloatType const& max_d_sq_upper_bound,
        af::tiny<FloatType, 3> const& grid_n_f,
        fractional<FloatType> const& coor_frac,
        gaussian_fourier_transformed<FloatType> const& gaussian_ft,
        int i_bv,
        int enlargement_step,
        int box_bound,
        int box_limit)
      {
        fractional<FloatType> d_frac(0,0,0);
        bool just_enlarged = false;
        while (true) {
          d_frac[i_bv] = enlargement_step
                       * (box_limit/grid_n_f[i_bv] - coor_frac[i_bv]);
          FloatType d_sq = unit_cell.length_sq(d_frac);
          FloatType abs_rho = scitbx::fn::absolute(gaussian_ft.rho_real(d_sq));
          if (abs_rho < rho_cutoff) {
            if (just_enlarged) {
              box_limit -= enlargement_step;
              break;
            }
            if (box_limit == box_bound) break;
            box_limit -= enlargement_step;
          }
          else {
            if (d_sq > max_d_sq_upper_bound) {
              excessive_radius = true;
              if (just_enlarged) box_limit -= enlargement_step;
              break;
            }
            scitbx::math::update_max(max_d_sq, d_sq);
            box_limit += enlargement_step;
            just_enlarged = true;
          }
        }
        return box_limit;
      }
    };

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
      typedef typename af::c_grid_padded_periodic<3> accessor_type;
      typedef typename accessor_type::index_type grid_point_type;
      typedef typename grid_point_type::value_type grid_point_element_type;
      typedef std::vector<FloatType> u_radius_cache_t;
      typedef std::vector<scitbx::sym_mat3<FloatType> > u_cart_cache_t;

      typedef af::versa<FloatType, accessor_type>
        real_map_type;
      typedef af::versa<std::complex<FloatType>, accessor_type>
        complex_map_type;

      sampling_base() {}

      sampling_base(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<XrayScattererType> const& scatterers,
        xray::scattering_type_registry const& scattering_type_registry,
        FloatType const& u_base,
        FloatType const& wing_cutoff,
        FloatType const& exp_table_one_over_step_size,
        FloatType const& tolerance_positive_definite,
        bool use_u_base_as_u_extra=false);

      uctbx::unit_cell const&
      unit_cell() { return unit_cell_; }

      FloatType
      u_base() const { return u_base_; }

      FloatType
      u_extra() const { return u_extra_; }

      FloatType
      u_min() const { return u_min_; }

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

      std::size_t
      max_sampling_box_n_points() const { return max_sampling_box_n_points_; }

      std::size_t
      sum_sampling_box_n_points() const { return sum_sampling_box_n_points_; }

      double
      ave_sampling_box_n_points() const
      {
        if (n_contributing_scatterers_ == 0) return 0;
        return static_cast<double>(sum_sampling_box_n_points_)
             / static_cast<double>(n_contributing_scatterers_);
      }

      grid_point_type const&
      max_sampling_box_edges() const { return max_sampling_box_edges_; }

      fractional<FloatType>
      max_sampling_box_edges_frac() const
      {
        CCTBX_ASSERT(map_accessor_.focus_size_1d() != 0);
        fractional<FloatType> r;
        for(std::size_t i=0;i<3;i++) {
          r[i] = FloatType(max_sampling_box_edges_[i])
               / map_accessor_.focus()[i];
        }
        return r;
      }

      af::shared<std::size_t> const&
      excessive_sampling_radius_i_seqs() const
      {
        return excessive_sampling_radius_i_seqs_;
      }

    protected:
      uctbx::unit_cell unit_cell_;
      std::size_t n_scatterers_passed_;
      FloatType u_base_;
      FloatType wing_cutoff_;
      FloatType exp_table_one_over_step_size_;
      FloatType tolerance_positive_definite_;
      std::size_t n_contributing_scatterers_;
      std::size_t n_anomalous_scatterers_;
      bool anomalous_flag_;
      FloatType u_min_;
      FloatType u_extra_;
      FloatType rho_cutoff_;
      FloatType max_d_sq_upper_bound_;
      std::size_t exp_table_size_;
      std::size_t max_sampling_box_n_points_;
      std::size_t sum_sampling_box_n_points_;
      grid_point_type max_sampling_box_edges_;
      af::shared<std::size_t> excessive_sampling_radius_i_seqs_;
      accessor_type map_accessor_;
      u_radius_cache_t u_radius_cache_;
      u_cart_cache_t u_cart_cache_;

      void
      update_sampling_box_statistics(
        std::size_t n_points,
        grid_point_type const& box_edges)
      {
        scitbx::math::update_max(max_sampling_box_n_points_, n_points);
        sum_sampling_box_n_points_ += n_points;
        detail::array_update_max(max_sampling_box_edges_, box_edges);
      }
  };

  template <typename FloatType,
            typename XrayScattererType>
  sampling_base<FloatType, XrayScattererType>
  ::sampling_base(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<XrayScattererType> const& scatterers,
    xray::scattering_type_registry const& scattering_type_registry,
    FloatType const& u_base,
    FloatType const& wing_cutoff,
    FloatType const& exp_table_one_over_step_size,
    FloatType const& tolerance_positive_definite,
    bool use_u_base_as_u_extra)
  :
    unit_cell_(unit_cell),
    n_scatterers_passed_(scatterers.size()),
    u_base_(u_base),
    wing_cutoff_(wing_cutoff),
    exp_table_one_over_step_size_(exp_table_one_over_step_size),
    tolerance_positive_definite_(tolerance_positive_definite),
    n_contributing_scatterers_(0),
    n_anomalous_scatterers_(0),
    anomalous_flag_(false),
    u_min_(-1),
    u_extra_(-1),
    rho_cutoff_(1),
    max_d_sq_upper_bound_(unit_cell.shortest_vector_sq()),
    exp_table_size_(0),
    max_sampling_box_n_points_(0),
    sum_sampling_box_n_points_(0),
    max_sampling_box_edges_(0,0,0)
  {
    CCTBX_ASSERT(u_base >= 0);
    scitbx::mat3<FloatType> orth_mx = unit_cell_.orthogonalization_matrix();
    if (orth_mx[3] != 0 || orth_mx[6] != 0 || orth_mx[7] != 0) {
      throw error(
        "Fatal Programming Error:"
        " Real-space sampling of model electron density"
        " is optimized for orthogonalization matrix"
        " according to the PDB convention. The orthogonalization"
        " matrix passed is not compatible with this convention.");
    }
    u_radius_cache_.reserve(scatterers.size());
    u_cart_cache_.reserve(scatterers.size());
    bool have_u_min = false;
    for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
      XrayScattererType const& scatterer = scatterers[i_seq];
      FloatType w = scatterer.weight();
      if (w == 0) continue;
      n_contributing_scatterers_++;
      if (scatterer.fdp != 0) {
        n_anomalous_scatterers_++;
      }
      FloatType u_min=0.;
      if (scatterer.flags.use_u_aniso()) {
        scitbx::sym_mat3<FloatType>
            u_cart = adptbx::u_star_as_u_cart(unit_cell_, scatterer.u_star);
        if (scatterer.flags.use_u_iso()) {
          u_cart[0] += scatterer.u_iso;
          u_cart[1] += scatterer.u_iso;
          u_cart[2] += scatterer.u_iso;
        }
        u_cart_cache_.push_back(u_cart);
        scitbx::vec3<FloatType>
          u_cart_eigenvalues = adptbx::eigenvalues(u_cart_cache_.back());
        u_radius_cache_.push_back(af::max(u_cart_eigenvalues));
        u_min = af::min(u_cart_eigenvalues);
      }
      else if (scatterer.flags.use_u_iso()) {
        u_radius_cache_.push_back(scatterer.u_iso);
        u_min = scatterer.u_iso;
      }
      else {
        u_radius_cache_.push_back(0);
        u_min = 0;
      }
      if (!have_u_min) {
        u_min_ = u_min;
        have_u_min = true;
      }
      else {
        u_min_ = std::min(u_min_, u_min);
      }
    }
    CCTBX_ASSERT(n_contributing_scatterers_ == u_radius_cache_.size());
    if (have_u_min) {
      if (!use_u_base_as_u_extra) {
        u_extra_ = u_base_ - u_min_;
      }
      else {
        u_extra_ = u_base_;
        u_base_ += u_min_;
      }
    }
    // determine rho_cutoff as average rho(0) * wing_cutoff
    {
      // exp function is not actually called by gaussian_ft.rho_real_0() below.
      detail::exponent_table<FloatType> exp_table(0);
      FloatType sum_abs_rho_real_0 = 0;
      typename u_radius_cache_t::const_iterator u_radius =
        u_radius_cache_.begin();
      typename u_cart_cache_t::const_iterator u_cart =
        u_cart_cache_.begin();
      detail::gaussian_fourier_transformed<FloatType> gaussian_ft;
      for(std::size_t i_seq=0;i_seq<scatterers.size();i_seq++) {
        XrayScattererType const& scatterer = scatterers[i_seq];
        FloatType w = scatterer.weight();
        if (w == 0) continue;
        eltbx::xray_scattering::gaussian const&
          gaussian = scattering_type_registry.gaussian_not_optional(
            scatterer.scattering_type);
        if (scatterer.flags.use_u_aniso()) {
          u_radius++;
          gaussian_ft = detail::gaussian_fourier_transformed<FloatType>(
            exp_table,
            gaussian, scatterer.fp, scatterer.fdp, w,
            *u_cart++, u_extra_);
        }
        else {
          gaussian_ft = detail::gaussian_fourier_transformed<FloatType>(
            exp_table,
            gaussian, scatterer.fp, scatterer.fdp, w,
            *u_radius++, u_extra_);
        }
        sum_abs_rho_real_0 += scitbx::fn::absolute(gaussian_ft.rho_real_0());
      }
      CCTBX_ASSERT(u_radius == u_radius_cache_.end());
      CCTBX_ASSERT(u_cart == u_cart_cache_.end());
      if (sum_abs_rho_real_0 != 0) {
        rho_cutoff_ = sum_abs_rho_real_0 / n_contributing_scatterers_
                    * wing_cutoff_;
      }
    }
    CCTBX_ASSERT(rho_cutoff_ > 0);
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SAMPLING_BASE_H
