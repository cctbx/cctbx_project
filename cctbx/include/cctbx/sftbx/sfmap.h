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
#include <cctbx/array_family/small.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/reductions.h>
#include <cctbx/maps/accessors.h>
#include <cctbx/maps/gridding.h>
#include <cctbx/maps/sym_tags.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_asu.h>

#include <cctbx/array_family/simple_io.h> // XXX

#if (defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
#include <cctbx/sftbx/xray_scatterer.h>
#endif

namespace cctbx { namespace sftbx {

  namespace sfmap {

    typedef maps::padded_grid_p1<3> accessor_type;
    typedef accessor_type::index_type grid_point_type;
    typedef grid_point_type::value_type grid_point_element_type;

  } // namespace detail

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

    template <typename FloatType>
    struct calc_shell
    {
      sfmap::grid_point_type radii;
      FloatType max_d_sq;

      calc_shell(
        const uctbx::UnitCell& ucell,
        const FloatType& wing_cutoff,
        const sfmap::grid_point_type& grid_n,
        const af::small<FloatType, 6>& ae,
        const af::small<FloatType, 6>& be,
        exponent_table<FloatType>& exp_table)
        : max_d_sq(0)
      {
        af::tiny<FloatType, 3> grid_n_f = grid_n;
        FloatType rho_cutoff = af::sum(ae.const_ref()) * wing_cutoff;
        for(std::size_t i_basis_vec=0;i_basis_vec<3;i_basis_vec++) {
          sfmap::grid_point_element_type ig;
          for (ig = 1; ig < grid_n[i_basis_vec]; ig++) {
            fractional<FloatType> d_frac(0,0,0);
            d_frac[i_basis_vec] = ig / grid_n_f[i_basis_vec];
            FloatType d_sq = ucell.Length2(d_frac);
            FloatType rho_g = 0;
            for (std::size_t i=0;i<ae.size();i++) {
              rho_g += ae[i] * exp_table(be[i] * d_sq);
            }
            if (rho_g < rho_cutoff) {
              break;
            }
            if (max_d_sq < d_sq) max_d_sq = d_sq;
          }
          if (ig == grid_n[i_basis_vec]) {
            throw error(
              "Excessive radius for real-space sampling"
              " of model electron density.");
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
      typedef af::versa<FloatType, sfmap::accessor_type> map_type;

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
        const FloatType& exp_table_one_over_step_size = -100)
        : ucell_(ucell),
          map_(sfmap::accessor_type(grid_logical, grid_physical)),
          u_extra_(u_extra),
          wing_cutoff_(wing_cutoff),
          exp_table_one_over_step_size_(exp_table_one_over_step_size),
          exp_table_size_(0),
          max_shell_radii_(0,0,0)
      {
        detail::exponent_table<FloatType> exp_table(
          exp_table_one_over_step_size);
        af::tiny<FloatType, 9> orth_mx = ucell_.getOrthogonalizationMatrix();
        if (orth_mx[3] != 0 || orth_mx[6] != 0 || orth_mx[7] != 0) {
          throw error(
            "Fatal Programming Error:"
            " Function is hand-optimized for orthogonalization matrix according"
            " to the PDB convention. The orthogonalization matrix passed is"
            " not compatible with this convention.");
        }
        for(
#if !((defined(BOOST_MSVC) && BOOST_MSVC <= 1300)) // VC++ 7.0
            const XrayScattererType*
#else
            const sftbx::XrayScatterer<double, eltbx::CAASF_WK1995>*
#endif
            site=sites.begin();site!=sites.end();site++)
        {
          using constants::four_pi;
          using constants::four_pi_sq;
          if (site->w() == 0) continue;
          cctbx_assert(site->fpfdp().imag() == 0);
          cctbx_assert(!site->isAnisotropic());
          const fractional<FloatType>& coor_frac = site->Coordinates();
          FloatType b_incl_extra = adptbx::U_as_B(site->Uiso() + u_extra_);
          // Calculate reciprocal space "form factors"
          af::small<FloatType, 6> ae;
          af::small<FloatType, 6> be;
          std::size_t i;
          for(i=0;i<site->CAASF().n_ab();i++) {
            FloatType b_i_b_incl_extra = site->CAASF().b(i) + b_incl_extra;
            FloatType f = std::pow(four_pi / b_i_b_incl_extra, 1.5);
            ae.push_back(site->w() * site->CAASF().a(i) * f);
            be.push_back(-four_pi_sq / b_i_b_incl_extra);
          }
          FloatType f = std::pow(four_pi / b_incl_extra, 1.5);
          ae.push_back(site->w() * site->CAASF().c() * f);
          be.push_back(-four_pi_sq / b_incl_extra);
          const sfmap::grid_point_type& grid_nl = map_.accessor().n_logical();
          const sfmap::grid_point_type& grid_np = map_.accessor().n_physical();
          detail::calc_shell<FloatType> shell(
            ucell_, wing_cutoff_, grid_nl, ae, be, exp_table);
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
          FloatType* map_begin = map_.begin();
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
            FloatType rho_g(0);
            for (std::size_t i=0;i<ae.size();i++) {
              rho_g += ae[i] * exp_table(be[i] * d_sq);
            }
            map_begin[g0112 + detail::mod_positive(gp[2],grid_nl[2])] += rho_g;
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
      std::size_t exp_table_size() const { return exp_table_size_; }
      const sfmap::grid_point_type& max_shell_radii() const {
        return max_shell_radii_;
      }
      const map_type& map() const { return map_; }
            map_type& map()       { return map_; }

      template <typename TagType>
      void
      apply_symmetry(const maps::grid_tags<TagType>& tags)
      {
        tags.sum_sym_equiv_points(map_.ref());
      }

      void
      eliminate_u_extra_and_normalize(
        const af::const_ref<Miller::Index>& miller_indices,
        af::ref<std::complex<FloatType> > structure_factors) const
      {
        FloatType norm = ucell_.getVolume()
                       / af::product(map_.accessor().n_logical().const_ref());
        eliminate_u_extra(
          ucell_, u_extra_, miller_indices, structure_factors, norm);
      }

    private:
      uctbx::UnitCell ucell_;
      map_type map_;
      FloatType u_extra_;
      FloatType wing_cutoff_;
      FloatType exp_table_one_over_step_size_;
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
    const af::const_ref<FloatType>& transformed_real_map,
    const IndexType& n_complex,
    bool conjugate = true)
  {
    cctbx_assert(
         transformed_real_map.size()
      >= 2 * af::compile_time_product<3>::get(n_complex));
    af::const_ref<std::complex<FloatType> > complex_map(
      reinterpret_cast<const std::complex<FloatType>*>(
        transformed_real_map.begin()), transformed_real_map.size()/2);
    sgtbx::ReciprocalSpaceASU asu(sginfo);
    af::shared<Miller::Index> miller_indices;
    af::shared<std::complex<FloatType> > structure_factors;
    Miller::Index h;
    Miller::Index mh;
    IndexType loop_i;
    std::size_t map_i = 0;
    for(loop_i[0] = 0; loop_i[0] < n_complex[0]; loop_i[0]++) {
      h[0] = maps::ih_as_h(loop_i[0], n_complex[0]);
      mh[0] = -h[0];
    for(loop_i[1] = 0; loop_i[1] < n_complex[1]; loop_i[1]++) {
      h[1] = maps::ih_as_h(loop_i[1], n_complex[1]);
      mh[1] = -h[1];
    for(loop_i[2] = 0; loop_i[2] < n_complex[2]; loop_i[2]++, map_i++) {
      h[2] = loop_i[2];
      mh[2] = -h[2];
      if (ucell.Q(h) <= max_q && map_i) {
        cctbx_assert(loop_i[0] != n_complex[0]/2);
        cctbx_assert(loop_i[1] != n_complex[1]/2);
        if (asu.isInASU(h)) {
          if (!sginfo.SgOps().isSysAbsent(h)) {
            miller_indices.push_back(h);
            if (!conjugate) {
              structure_factors.push_back(complex_map[map_i]);
            }
            else {
              structure_factors.push_back(std::conj(complex_map[map_i]));
            }
          }
        }
        else if (h[2] && asu.isInASU(mh)) {
          if (!sginfo.SgOps().isSysAbsent(h)) {
            miller_indices.push_back(mh);
            if (!conjugate) {
              structure_factors.push_back(std::conj(complex_map[map_i]));
            }
            else {
              structure_factors.push_back(complex_map[map_i]);
            }
          }
        }
      }
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

#endif // CCTBX_XRAY_SCATTERER_H
