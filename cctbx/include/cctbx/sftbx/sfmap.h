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

  namespace detail {

    template <typename NumType>
    inline
    NumType
    pow2(const NumType& x) { return x * x; }

    template <typename FloatType>
    inline
    FloatType
    add_for_rounding(const FloatType& x) {
      if (x < FloatType(0)) return x - FloatType(.5);
      return x + FloatType(.5);
    }

    template <typename FloatType,
              typename IndexType>
    IndexType
    calc_nearest_grid_point(const fractional<FloatType>& coor,
                            const IndexType& grid_n)
    {
      typedef typename IndexType::value_type index_value_type;
      IndexType grid_point;
      for(std::size_t i=0;i<3;i++) {
        grid_point[i] = index_value_type(
          add_for_rounding(coor[i] * grid_n[i]));
      }
      return grid_point;
    }

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
        cctbx_assert(xs >= 0); // Use NDEBUG
        std::size_t i(xs + FloatType(.5));
        if (i >= table_.size()) expand(i + 1);
        return table_[i];
      }
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
  class sampled_model_density
  {
    public:
      typedef af::versa<FloatType, maps::padded_grid_p1<3> > map_type;
      typedef typename map_type::accessor_type map_accessor_type;
      typedef typename map_accessor_type::index_type grid_point_type;

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
        const FloatType& max_q,
        const FloatType& grid_resolution_factor,
        const FloatType& u_extra,
        const grid_point_type& grid_logical,
        const grid_point_type& grid_physical,
        const FloatType& wing_cutoff = 0.01,
        const FloatType& exp_table_one_over_step_size = -1000)
        : ucell_(ucell),
          u_extra_(u_extra),
          map_(map_accessor_type(grid_logical, grid_physical)),
          max_q_(max_q),
          grid_resolution_factor_(grid_resolution_factor),
          wing_cutoff_(wing_cutoff),
          exp_table_(exp_table_one_over_step_size)
      {
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
          cctbx_assert(site->fpfdp().imag() == 0);
          cctbx_assert(!site->isAnisotropic());
          const fractional<FloatType>& coor_frac = site->Coordinates();
          cartesian<FloatType> coor_cart = ucell.orthogonalize(coor_frac);
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
          FloatType max_d_sq = calc_max_sampling_distance_sq(ae, be);
          //XXXstd::cout << "max_d_sq: " << max_d_sq << std::endl;
          //Calculate limits of shell search
          cartesian<FloatType> max_d_cart;
          max_d_cart.fill(std::sqrt(max_d_sq));
          //XXXstd::cout << "max_d_cart: " << max_d_cart.const_ref() << std::endl;
          fractional<FloatType> max_d_frac = ucell.fractionalize(max_d_cart);
          //XXXstd::cout << "max_d_frac: " << max_d_frac.const_ref() << std::endl;
          const grid_point_type& grid_logical = map_.accessor().n_logical();
          //XXXstd::cout << "grid_logical: " << grid_logical.const_ref() << std::endl;
          af::int3 shell_limit;
          for(i=0;i<3;i++) {
            //Round number down to nearest integer as you will never "make it"
            //to the grid + 1 point
            shell_limit[i] = int(
              std::floor(max_d_frac[i] * grid_logical[i]) + .5);
          }
          //XXXstd::cout << "shell_limit: " << shell_limit.const_ref() << std::endl;
          grid_point_type pivot = detail::calc_nearest_grid_point(
            coor_frac, grid_logical);
          //XXXstd::cout << "pivot: " << pivot.const_ref() << std::endl;
          // XXX TODO pre-compute pivot - coor_frac
          grid_point_type ip;
          fractional<FloatType> g_frac;
          for(ip[0] = -shell_limit[0]; ip[0] <= shell_limit[0]; ip[0]++) {
            g_frac[0] = FloatType(pivot[0] + ip[0]) / grid_logical[0];
          for(ip[1] = -shell_limit[1]; ip[1] <= shell_limit[1]; ip[1]++) {
            g_frac[1] = FloatType(pivot[1] + ip[1]) / grid_logical[1];
          for(ip[2] = -shell_limit[2]; ip[2] <= shell_limit[2]; ip[2]++) {
            g_frac[2] = FloatType(pivot[2] + ip[2]) / grid_logical[2];
            fractional<FloatType> d_frac = g_frac - coor_frac;
            FloatType d_sq = ucell.Length2(d_frac);
            //XXXstd::cout << "ip: " << ip.const_ref() << std::endl;
            //XXXstd::cout << "d_frac: " << d_frac.const_ref() << std::endl;
            //XXXstd::cout << "d_sq: " << d_sq << std::endl;
            if (d_sq > max_d_sq) continue;
            FloatType rho_g(0);
            for (std::size_t i=0;i<ae.size();i++) {
              rho_g += ae[i] * exp_table_(be[i] * d_sq);
            }
            //std::cout << "rho_g: " << rho_g << std::endl;
            map_(pivot + ip) += rho_g;
          }}}
        }
      }

      const uctbx::UnitCell& ucell() { return ucell_; }
      const FloatType& max_q() const { return max_q_; }
      const FloatType& grid_resolution_factor() const {
        return grid_resolution_factor_;
      }
      const FloatType& u_extra() const { return u_extra_; }
      const FloatType& wing_cutoff() const { return wing_cutoff_; }
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
      FloatType max_q_;
      FloatType grid_resolution_factor_;
      FloatType wing_cutoff_;
      FloatType u_extra_;
      exponent_table<FloatType> exp_table_;

      FloatType
      calc_max_sampling_distance_sq(
        const af::small<FloatType, 6>& ae,
        const af::small<FloatType, 6>& be)
      {
        FloatType rho_origin = af::sum(ae.const_ref());
        FloatType sampling_unit = grid_resolution_factor_ / std::sqrt(max_q_);
        for (std::size_t radius_step = 1;; radius_step++) {
          FloatType d_sq = detail::pow2(radius_step * sampling_unit);
          FloatType rho_g = 0;
          for (std::size_t i=0;i<ae.size();i++) {
            rho_g += ae[i] * exp_table_(be[i] * d_sq);
          }
          if (rho_g < rho_origin * wing_cutoff_) {
            return d_sq;
          }
          cctbx_assert(
            radius_step < af::max(map_.accessor().n_logical().const_ref()));
        }
      }
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

  typedef af::grid<3>::index_type grid_point_type;
  typedef grid_point_type::value_type grid_point_element_type;

  template <typename nType>
  grid_point_type
  h_as_ih_array(const Miller::Index& h, const nType& n)
  {
    grid_point_type ih;
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
        grid_point_type ih = h_as_ih_array(h, n_complex);
        cctbx_assert(ih >= grid_point_element_type(0));
        if (!conjugate) {
          map(ih) = semi.ShiftPhase(e, structure_factors(i));
        }
        else {
          map(ih) = std::conj(semi.ShiftPhase(e, structure_factors(i)));
        }
        if (h[2] == 0 && !semi.isCentric()) {
          grid_point_type ih = h_as_ih_array(-h, n_complex);
          cctbx_assert(ih >= grid_point_element_type(0));
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
