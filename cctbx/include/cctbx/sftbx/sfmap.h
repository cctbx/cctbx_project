// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Copyright (c) 2002 Airlie McCoy.

   Revision history:
     2002 May: Created based on phaser/src/MapFFT.cc by Airlie McCoy (rwgk)
 */

#ifndef CCTBX_SFMAP_H
#define CCTBX_SFMAP_H

#include <vector>
#include <cctbx/error.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/reductions.h>
#include <cctbx/maps/accessors.h>
#include <cctbx/maps/gridding.h>
#include <cctbx/maps/sym_tags.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_asu.h>

#if (defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
#include <cctbx/sftbx/xray_scatterer.h>
#endif

namespace cctbx { namespace sftbx {

  namespace sfmap {

    typedef maps::padded_grid_p1<3> accessor_type;
    typedef accessor_type::index_type grid_point_type;
    typedef grid_point_type::value_type grid_point_element_type;

  }

  //! Artificial temperature factor for the treatment of aliasing problems.
  /*! Reference:

        Gerard Bricogne (2001),
        International Tables for Crystallography, Volume B, 2001, p. 87
        (end of section 1.3.4.4.5).

      @param max_q = 1/(d*max)^2
      @param grid_resolution_factor = 1/(2*sigma)
      @param quality_factor = Q
      @param max_u_extra is a user-defined upper limit.

      max_q can be obtained with cctbx::uctbx::UnitCell::Q().

      quality_factor = 100 for 1% accuracy.
    */
  template <typename FloatType>
  FloatType
  calc_u_extra(
    const FloatType& max_q,
    const FloatType& grid_resolution_factor,
    const FloatType& quality_factor = 100,
    const FloatType& max_u_extra = adptbx::B_as_U(1000))
  {
    FloatType numerator = adptbx::B_as_U(std::log10(quality_factor));
    FloatType sigma = 1 / (2 * grid_resolution_factor);
    FloatType denominator = sigma * (sigma - 1) * max_q;
    if (max_u_extra * denominator > numerator) {
      return numerator / denominator;
    }
    return max_u_extra;
  }

  template <typename FloatType>
  void
  eliminate_u_extra(
    const uctbx::UnitCell& ucell,
    const FloatType& u_extra,
    const af::const_ref<Miller::Index>& miller_indices,
    af::ref<std::complex<FloatType> > structure_factors,
    const FloatType& norm = 1)
  {
    cctbx_assert(miller_indices.size() == structure_factors.size());
    //    f *= exp(b_extra * d_star_sq(h)/4)
    // => f *= exp(2*pi*pi*u_extra * d_star_sq(h))
    FloatType tppu = cctbx::constants::two_pi_sq * u_extra;
    for(std::size_t i=0;i<miller_indices.size();i++) {
      structure_factors[i] *= norm * FloatType(std::exp(
        tppu * ucell.Q(miller_indices[i])));
    }
  }

  namespace detail {

    template <typename NumType>
    inline
    NumType
    abs(const NumType& x) {
      if (x < NumType(0)) return -x;
      return x;
    }

    // self-expanding exponent table
    template <typename FloatType>
    class exponent_table
    {
      public:
        exponent_table() {}
        explicit
        exponent_table(const FloatType& one_over_step_size)
          : one_over_step_size_(one_over_step_size)
        {}
        FloatType operator()(const FloatType& x)
        {
          if (one_over_step_size_ == 0) return std::exp(x);
          FloatType xs = x * one_over_step_size_;
          cctbx_assert(xs >= 0);
          std::size_t i(xs + FloatType(.5));
          if (i >= table_.size()) expand(i + 1);
          return table_[i];
        }
        const std::vector<FloatType>& table() const { return table_; }
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

    template <typename FloatType,
              typename CaasfType>
    class caasf_fourier_transformed
    {
      public:
        caasf_fourier_transformed(
          const CaasfType& caasf,
          const std::complex<FloatType>& fpfdp,
          const FloatType& w,
          const FloatType& u_incl_extra)
        {
          using constants::four_pi;
          using constants::four_pi_sq;
          FloatType b_incl_extra = adptbx::U_as_B(u_incl_extra);
          std::size_t i = 0;
          for(;i<caasf.n_ab();i++) {
            FloatType b_i_b_incl_extra = caasf.b(i) + b_incl_extra;
            FloatType f = std::pow(four_pi / b_i_b_incl_extra, 1.5);
            a_real_[i] = w * caasf.a(i) * f;
            b_real_[i] = -four_pi_sq / b_i_b_incl_extra;
          }
          FloatType f = std::pow(four_pi / b_incl_extra, 1.5);
          a_real_[i] = w * (caasf.c() + fpfdp.real()) * f;
          b_real_[i] = -four_pi_sq / b_incl_extra;
          a_imag_ = w * fpfdp.imag() * f;
          b_imag_ = b_real_[i];
        }

        FloatType rho_real_0() const
        {
          return af::sum(a_real_.const_ref());
        }

        FloatType
        rho_real(exponent_table<FloatType>& exp_table,
                 const FloatType& d_sq) const
        {
          FloatType r(0);
          for (std::size_t i=0;i<a_real_.size();i++) {
            r += a_real_[i] * exp_table(b_real_[i] * d_sq);
          }
          return r;
        }

        FloatType
        rho_imag(exponent_table<FloatType>& exp_table,
                 const FloatType& d_sq) const
        {
          return a_imag_ * exp_table(b_imag_ * d_sq);
        }

      private:
        af::tiny<FloatType, CaasfType::n_plus_1> a_real_;
        af::tiny<FloatType, CaasfType::n_plus_1> b_real_;
        FloatType a_imag_;
        FloatType b_imag_;
    };

    template <typename FloatType>
    struct calc_shell
    {
      sfmap::grid_point_type radii;
      FloatType max_d_sq;

      template <typename CaasfType>
      calc_shell(
        const uctbx::UnitCell& ucell,
        const FloatType& wing_cutoff,
        const sfmap::grid_point_type& grid_n,
        const caasf_fourier_transformed<FloatType, CaasfType>& caasf_ft,
        exponent_table<FloatType>& exp_table)
        : max_d_sq(0)
      {
        af::tiny<FloatType, 3> grid_n_f = grid_n;
        FloatType rho_cutoff = detail::abs(
          caasf_ft.rho_real_0() * wing_cutoff);
        for(std::size_t i_basis_vec=0;i_basis_vec<3;i_basis_vec++) {
          sfmap::grid_point_element_type ig;
          for (ig = 1; ig < grid_n[i_basis_vec]; ig++) {
            fractional<FloatType> d_frac(0,0,0);
            d_frac[i_basis_vec] = ig / grid_n_f[i_basis_vec];
            FloatType d_sq = ucell.Length2(d_frac);
            if (detail::abs(caasf_ft.rho_real(exp_table, d_sq)) < rho_cutoff) {
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
    add_for_rounding(const FloatType& x) {
      if (x < FloatType(0)) return x - FloatType(.5);
      return x + FloatType(.5);
    }

    template <typename FloatType>
    inline
    sfmap::grid_point_type
    calc_nearest_grid_point(const fractional<FloatType>& coor,
                            const sfmap::grid_point_type& grid_n)
    {
      sfmap::grid_point_type grid_point;
      for(std::size_t i=0;i<3;i++) {
        grid_point[i] = sfmap::grid_point_element_type(
          add_for_rounding(coor[i] * grid_n[i]));
      }
      return grid_point;
    }

    template <typename ArrayMaxType,
              typename ArrayValType>
    ArrayMaxType&
    array_update_max(ArrayMaxType& am, const ArrayValType& av)
    {
      for(std::size_t i=0;i<am.size();i++) {
        if (am[i] < av[i]) am[i] = av[i];
      }
      return am;
    }

    template <typename IntegerType1,
              typename IntegerType2>
    inline
    IntegerType1
    mod_positive(IntegerType1 ix, const IntegerType2& iy)
    {
      ix %= iy;
      if (ix < 0) ix += iy;
      return ix;
    }

  } // namespace detail

  template <typename FloatType>
  class sampled_model_density
  {
    public:
      typedef af::versa<FloatType, sfmap::accessor_type>
        map_real_type;
      typedef af::versa<std::complex<FloatType>, sfmap::accessor_type>
        map_complex_type;

      sampled_model_density() {}

#if !((defined(BOOST_MSVC) && BOOST_MSVC <= 1300)) // VC++ 7.0
      template <typename XrayScattererType>
#endif
      sampled_model_density(
        const uctbx::UnitCell& ucell,
#if !((defined(BOOST_MSVC) && BOOST_MSVC <= 1300)) // VC++ 7.0
        const af::shared<XrayScattererType>& sites,
#else
        const af::shared<
          sftbx::XrayScatterer<double, eltbx::CAASF_WK1995> >& sites,
#endif
        const sfmap::grid_point_type& grid_logical,
        const sfmap::grid_point_type& grid_physical,
        const FloatType& u_extra = 0.25,
        const FloatType& wing_cutoff = 1.e-3,
        const FloatType& exp_table_one_over_step_size = -100,
        bool force_complex = false,
        bool electron_density_must_be_positive = true)
        : ucell_(ucell),
          n_passed_scatterers_(sites.size()),
          u_extra_(u_extra),
          wing_cutoff_(wing_cutoff),
          exp_table_one_over_step_size_(exp_table_one_over_step_size),
          n_contributing_scatterers_(0),
          n_anomalous_scatterers_(0),
          exp_table_size_(0),
          max_shell_radii_(0,0,0)
      {
        typedef
#if !((defined(BOOST_MSVC) && BOOST_MSVC <= 1300)) // VC++ 7.0
          XrayScattererType
#else
          sftbx::XrayScatterer<double, eltbx::CAASF_WK1995>
#endif
            xray_scatterer_type;
        typedef typename
          xray_scatterer_type::caasf_type::base_type
            caasf_type;
        af::tiny<FloatType, 9> orth_mx = ucell_.getOrthogonalizationMatrix();
        if (orth_mx[3] != 0 || orth_mx[6] != 0 || orth_mx[7] != 0) {
          throw error(
            "Fatal Programming Error:"
            " Real-space sampling of model electron density"
            " is hand-optimized for orthogonalization matrix"
            " according to the PDB convention. The orthogonalization"
            " matrix passed is not compatible with this convention.");
        }
        const xray_scatterer_type* site;
        for(site=sites.begin();site!=sites.end();site++)
        {
          if (site->w() == 0) continue;
          n_contributing_scatterers_++;
          if (site->fpfdp().imag() != 0) {
            n_anomalous_scatterers_++;
          }
        }
        bool friedel_flag;
        FloatType* map_begin;
        if (n_anomalous_scatterers_ == 0 && !force_complex) {
          friedel_flag = true;
          map_accessor_ = sfmap::accessor_type(grid_logical, grid_physical);
          map_real_.resize(map_accessor_);
          map_begin = map_real_.begin();
        }
        else {
          friedel_flag = false;
          map_accessor_ = sfmap::accessor_type(grid_logical, grid_logical);
          map_complex_.resize(map_accessor_);
          map_begin = reinterpret_cast<FloatType*>(map_complex_.begin());
        }
        const sfmap::grid_point_type& grid_nl = map_accessor_.n_logical();
        const sfmap::grid_point_type& grid_np = map_accessor_.n_physical();
        detail::exponent_table<FloatType> exp_table(
          exp_table_one_over_step_size);
        for(site=sites.begin();site!=sites.end();site++)
        {
          if (site->w() == 0) continue;
          cctbx_assert(!site->isAnisotropic());
          FloatType fdp = site->fpfdp().imag();
          const fractional<FloatType>& coor_frac = site->Coordinates();
          detail::caasf_fourier_transformed<FloatType, caasf_type> caasf_ft(
            site->CAASF(), site->fpfdp(), site->w(), site->Uiso() + u_extra_);
          detail::calc_shell<FloatType> shell(
            ucell_, wing_cutoff_, grid_nl, caasf_ft, exp_table);
          if (electron_density_must_be_positive) {
            if (   caasf_ft.rho_real_0() < 0
                || caasf_ft.rho_real(exp_table, shell.max_d_sq) < 0) {
              throw error("Negative electron density at sampling point.");
            }
          }
          detail::array_update_max(max_shell_radii_, shell.radii);
          sfmap::grid_point_type pivot = detail::calc_nearest_grid_point(
            coor_frac, grid_nl);
          // highly hand-optimized loop over points in shell
          sfmap::grid_point_type g_min = pivot - shell.radii;
          sfmap::grid_point_type g_max = pivot + shell.radii;
          sfmap::grid_point_type gp;
          sfmap::grid_point_element_type g01, g0112;
          FloatType f0, f1, f2;
          FloatType c00, c01, c11, c02, c12, c22;
          for(gp[0] = g_min[0]; gp[0] <= g_max[0]; gp[0]++) {
            g01 = detail::mod_positive(gp[0],grid_nl[0]) * grid_np[1];
            f0 = FloatType(gp[0]) / grid_nl[0] - coor_frac[0];
            c00 = orth_mx[0] * f0;
          for(gp[1] = g_min[1]; gp[1] <= g_max[1]; gp[1]++) {
            g0112 = (g01+detail::mod_positive(gp[1],grid_nl[1])) * grid_np[2];
            f1 = FloatType(gp[1]) / grid_nl[1] - coor_frac[1];
            c01 = orth_mx[1] * f1 + c00;
            c11 = orth_mx[4] * f1;
          for(gp[2] = g_min[2]; gp[2] <= g_max[2]; gp[2]++) {
            f2 = FloatType(gp[2]) / grid_nl[2] - coor_frac[2];
            c02 = orth_mx[2] * f2 + c01;
            c12 = orth_mx[5] * f2 + c11;
            c22 = orth_mx[8] * f2;
            FloatType d_sq = c02*c02 + c12*c12 + c22*c22;
            if (d_sq > shell.max_d_sq) continue;
            std::size_t i_map = g0112 + detail::mod_positive(gp[2],grid_nl[2]);
            if (friedel_flag) {
              map_begin[i_map] += caasf_ft.rho_real(exp_table, d_sq);
            }
            else {
              i_map *= 2;
              map_begin[i_map] += caasf_ft.rho_real(exp_table, d_sq);
              if (fdp) {
                map_begin[i_map+1] += caasf_ft.rho_imag(exp_table, d_sq);
              }
            }
          }}}
        }
        exp_table_size_ = exp_table.table().size();
      }

      const uctbx::UnitCell& ucell() { return ucell_; }
      const FloatType& u_extra() const { return u_extra_; }
      const FloatType& wing_cutoff() const { return wing_cutoff_; }
      const FloatType& exp_table_one_over_step_size() const {
        return exp_table_one_over_step_size_;
      }
      std::size_t n_passed_scatterers() const {
        return n_passed_scatterers_;
      }
      std::size_t n_contributing_scatterers() const {
        return n_contributing_scatterers_;
      }
      std::size_t n_anomalous_scatterers() const {
        return n_anomalous_scatterers_;
      }
      bool friedel_flag() const { return map_real_.size() != 0; }
      const map_real_type& map_real() const { return map_real_; }
            map_real_type  map_real()       { return map_real_; }
      const map_complex_type& map_complex() const { return map_complex_; }
            map_complex_type  map_complex()       { return map_complex_; }
      std::size_t exp_table_size() const { return exp_table_size_; }
      const sfmap::grid_point_type& max_shell_radii() const {
        return max_shell_radii_;
      }
      fractional<FloatType>
      max_shell_radii_frac() const {
        fractional<FloatType> r;
        for(std::size_t i=0;i<3;i++) {
          r[i] = FloatType(max_shell_radii_[i]) / map_accessor_.n_logical()[i];
        }
        return r;
      }

      template <typename TagType>
      void
      apply_symmetry(const maps::grid_tags<TagType>& tags)
      {
        if (map_real_.size()) {
          tags.sum_sym_equiv_points(map_real_.ref());
        }
        else {
          tags.sum_sym_equiv_points(map_complex_.ref());
        }
      }

      void
      eliminate_u_extra_and_normalize(
        const af::const_ref<Miller::Index>& miller_indices,
        af::ref<std::complex<FloatType> > structure_factors) const
      {
        FloatType norm = ucell_.getVolume()
          / af::product(map_accessor_.n_logical().const_ref());
        eliminate_u_extra(
          ucell_, u_extra_, miller_indices, structure_factors, norm);
      }

    private:
      uctbx::UnitCell ucell_;
      std::size_t n_passed_scatterers_;
      FloatType u_extra_;
      FloatType wing_cutoff_;
      FloatType exp_table_one_over_step_size_;
      std::size_t n_anomalous_scatterers_;
      std::size_t n_contributing_scatterers_;
      sfmap::accessor_type map_accessor_;
      map_real_type map_real_;
      map_complex_type map_complex_;
      std::size_t exp_table_size_;
      sfmap::grid_point_type max_shell_radii_;
  };

  template <typename FloatType,
            typename IndexType>
  std::pair<
    af::shared<Miller::Index>,
    af::shared<std::complex<FloatType> > >
  collect_structure_factors(
    const uctbx::UnitCell& ucell,
    const sgtbx::SpaceGroupInfo& sginfo,
    const FloatType& max_q,
    const af::const_ref<std::complex<FloatType> >& complex_map,
    const IndexType& n_complex,
    bool friedel_flag,
    bool conjugate)
  {
    cctbx_assert(complex_map.size() >= af::product(n_complex.const_ref()));
    sgtbx::ReciprocalSpaceASU asu(sginfo);
    const sgtbx::SpaceGroup& sgops = sginfo.SgOps();
    af::shared<Miller::Index> miller_indices;
    af::shared<std::complex<FloatType> > structure_factors;
    Miller::Index h;
    Miller::Index mh;
    uctbx::incremental_d_star_sq<FloatType> incr_d_star_sq(ucell);
    IndexType loop_i;
    std::size_t map_i = 0;
    for(loop_i[0] = 0; loop_i[0] < n_complex[0]; loop_i[0]++) {
      h[0] = maps::ih_as_h(loop_i[0], n_complex[0]);
      mh[0] = -h[0];
      incr_d_star_sq.update0(h[0]);
    for(loop_i[1] = 0; loop_i[1] < n_complex[1]; loop_i[1]++) {
      h[1] = maps::ih_as_h(loop_i[1], n_complex[1]);
      mh[1] = -h[1];
      incr_d_star_sq.update1(h[1]);
    for(loop_i[2] = 0; loop_i[2] < n_complex[2]; loop_i[2]++, map_i++) {
      if (!friedel_flag) {
        h[2] = maps::ih_as_h(loop_i[2], n_complex[2]);
      }
      else {
        h[2] = loop_i[2];
      }
      if (incr_d_star_sq.get(h[2]) > max_q || map_i == 0) {
        continue;
      }
      cctbx_assert(loop_i[0] != n_complex[0]/2);
      cctbx_assert(loop_i[1] != n_complex[1]/2);
      if (!friedel_flag) {
        cctbx_assert(loop_i[2] != n_complex[2]/2);
      }
      mh[2] = -h[2];
      int asu_sign = asu.asu_sign(h, mh);
      if (asu_sign == 0) continue;
      sgtbx::sys_absent_test sa_test(sgops, h);
      if (sa_test.is_sys_absent()) continue;
      bool f_conj = false;
      if (friedel_flag) {
        if (asu_sign > 0) {
          miller_indices.push_back(h);
          f_conj = conjugate;
        }
        else {
          if (h[2] == 0) continue;
          miller_indices.push_back(mh);
          f_conj = !conjugate;
        }
      }
      else {
        if (((asu_sign < 0) != conjugate) && sa_test.is_centric()) continue;
        if (conjugate) miller_indices.push_back(mh);
        else           miller_indices.push_back(h);
      }
      if (f_conj) structure_factors.push_back(std::conj(complex_map[map_i]));
      else        structure_factors.push_back(complex_map[map_i]);
    }}}
    return std::make_pair(miller_indices, structure_factors);
  }

  template <typename nType>
  sfmap::grid_point_type
  h_as_ih_array(const Miller::Index& h, const nType& n)
  {
    sfmap::grid_point_type ih;
    const bool positive_only[] = {false, false, true};
    for(std::size_t i=0;i<3;i++) {
      ih[i] = maps::h_as_ih(h[i], n[i], positive_only[i]);
    }
    return ih;
  }

  // XXX needs to be updated to cope with anomalous data
  template <typename FloatType,
            typename IndexType>
  af::versa<std::complex<FloatType>, af::grid<3> >
  structure_factor_map(
    const sgtbx::SpaceGroup& sgops,
    const af::const_ref<Miller::Index>& miller_indices,
    const af::const_ref<std::complex<FloatType> >& structure_factors,
    const IndexType& n_complex,
    bool conjugate = true)
  {
    af::versa<std::complex<FloatType>, af::grid<3> > map(n_complex);
    for(std::size_t i=0;i<miller_indices.size();i++) {
      sgtbx::SymEquivMillerIndices semi = sgops.getEquivMillerIndices(
        miller_indices[i]);
      for(int e=0;e<semi.M(true);e++) {
        Miller::Index h = semi(e);
        if (h[2] < 0) continue;
        sfmap::grid_point_type ih = h_as_ih_array(h, n_complex);
        cctbx_assert(ih >= sfmap::grid_point_element_type(0));
        if (!conjugate) {
          map(ih) = semi.ShiftPhase(e, structure_factors(i));
        }
        else {
          map(ih) = std::conj(semi.ShiftPhase(e, structure_factors(i)));
        }
        if (h[2] == 0 && !semi.isCentric()) {
          sfmap::grid_point_type ih = h_as_ih_array(-h, n_complex);
          cctbx_assert(ih >= sfmap::grid_point_element_type(0));
          if (!conjugate) {
            map(ih) = std::conj(semi.ShiftPhase(e, structure_factors(i)));
          }
          else {
            map(ih) = semi.ShiftPhase(e, structure_factors(i));
          }
        }
      }
    }
    return map;
  }

}} // namespace cctbx::sftbx

#endif // CCTBX_SFMAP_H
