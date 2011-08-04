#ifndef CCTBX_MAPTBX_EIGHT_POINT_INTERPOLATION_H
#define CCTBX_MAPTBX_EIGHT_POINT_INTERPOLATION_H

#include <cctbx/coordinates.h>
#include <scitbx/math/modulo.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/interpolation.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace maptbx {

  template <typename IndexType,
            typename FloatType,
            typename SignedIntType=long>
  class get_corner
  {
    public:
      typedef typename IndexType::value_type iv_t;

      get_corner() {}

      get_corner(
        IndexType const& grid_n,
        fractional<FloatType> const& x_frac)
      {
        for(std::size_t i=0;i<3;i++) {
          FloatType xn = x_frac[i] * static_cast<FloatType>(grid_n[i]);
          SignedIntType ixn = scitbx::math::float_int_conversions<
            FloatType, SignedIntType>::ifloor(xn);
          i_grid[i] = scitbx::math::mod_positive(
            ixn, static_cast<SignedIntType>(grid_n[i]));
          weights_[i][1] = xn - static_cast<FloatType>(ixn);
          weights_[i][0] = 1 - weights_[i][1];
        }
      }

      get_corner(
        scitbx::mat3<FloatType> const& gridding_matrix,
        scitbx::vec3<FloatType> const& site_cart)
      {
        scitbx::vec3<FloatType> grid_float = gridding_matrix * site_cart;
        for(std::size_t i=0;i<3;i++) {
          SignedIntType ixn = scitbx::math::float_int_conversions<
            FloatType, SignedIntType>::ifloor(grid_float[i]);
          i_grid[i] = ixn;
          weights_[i][1] = grid_float[i] - static_cast<FloatType>(ixn);
          weights_[i][0] = 1 - weights_[i][1];
        }
      }

      get_corner(
        crystal::direct_space_asu::asu_mappings<FloatType> & am,
        IndexType const& grid_n,
        fractional<FloatType> const& x_frac)
      {
        cartesian<FloatType> const & x_cart = am.process(x_frac).mappings().back()[0].mapped_site();
        fractional<FloatType> new_x_frac = am.unit_cell().fractionalize(x_cart);
        FloatType epsilon = scitbx::math::floating_point_epsilon<FloatType>::get() * 10;
        for ( std::size_t i=0; i<3; ++i )
          if ( std::abs(new_x_frac[i]) < epsilon )
            new_x_frac[i] = 0;
        for(std::size_t i=0;i<3;i++) {
          FloatType xn = new_x_frac[i] * static_cast<FloatType>(grid_n[i]);
          SignedIntType ixn = scitbx::math::float_int_conversions<
            FloatType, SignedIntType>::ifloor(xn);
          i_grid[i] = ixn;
          weights_[i][1] = xn - static_cast<FloatType>(ixn);
          weights_[i][0] = 1 - weights_[i][1];
        }
      }

      FloatType
      weight(iv_t s0, iv_t s1, iv_t s2) const
      {
        return weights_[0][s0] * weights_[1][s1] * weights_[2][s2];
      }

      IndexType
      closest_grid_point(IndexType const& grid_n) const
      {
        IndexType result = i_grid;
        for(std::size_t i=0;i<3;i++) {
          if (weights_[i][1] > weights_[i][0]) {
            result[i]++;
            if (result[i] == grid_n[i]) {
              result[i] = 0;
            }
          }
        }
        return result;
      }

      IndexType i_grid;

    protected:
      FloatType weights_[3][2];
  };

  template <
    typename MapFloatType,
    typename SiteFloatType>
  MapFloatType
  eight_point_interpolation(
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac)
  {
    typedef af::c_grid_padded<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    MapFloatType result = 0;
    for(iv_t s0=0;s0<2;s0++) { iv_t i0 = (corner.i_grid[0] + s0) % grid_n[0];
    for(iv_t s1=0;s1<2;s1++) { iv_t i1 = (corner.i_grid[1] + s1) % grid_n[1];
    for(iv_t s2=0;s2<2;s2++) { iv_t i2 = (corner.i_grid[2] + s2) % grid_n[2];
      result += map(i0,i1,i2) * corner.weight(s0,s1,s2);
    }}}
    return result;
  }

  template <
    typename MapFloatType,
    typename SiteFloatType>
  MapFloatType
  eight_point_interpolation(
    af::const_ref<MapFloatType, af::c_grid<3> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac)
  {
    typedef af::c_grid<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    MapFloatType result = 0;
    for(iv_t s0=0;s0<2;s0++) { iv_t i0 = (corner.i_grid[0] + s0) % grid_n[0];
    for(iv_t s1=0;s1<2;s1++) { iv_t i1 = (corner.i_grid[1] + s1) % grid_n[1];
    for(iv_t s2=0;s2<2;s2++) { iv_t i2 = (corner.i_grid[2] + s2) % grid_n[2];
      result += map(i0,i1,i2) * corner.weight(s0,s1,s2);
    }}}
    return result;
  }

  // derived from http://www.paulinternet.nl/?page=bicubic
  // also see: http://en.wikipedia.org/wiki/Tricubic_interpolation
  template <
    typename MapFloatType,
    typename SiteFloatType>
  MapFloatType
  tricubic_interpolation(
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac)
  {
    using scitbx::math::interpolate_at_point;
    using namespace std;
    typedef af::c_grid_padded<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    af::tiny<MapFloatType, 4> p(0.0);
    af::tiny<SiteFloatType, 3> xn;
    for (unsigned k = 0; k < 3; k++) {
      if (x_frac[k] < 0) {
        xn[k] = fmod((1.-x_frac[k])*static_cast<SiteFloatType>(grid_n[k]), 1.0);
      } else {
        xn[k] = fmod(x_frac[k]*static_cast<SiteFloatType>(grid_n[k]), 1.0);
      }
    }
    for (int i = -1; i < 3; i++) {
      iv_t u = (corner.i_grid[0] + i) % grid_n[0];
      af::tiny<MapFloatType, 4> pp(0.0);
      for (int j = -1; j < 3; j++) {
        iv_t v = (corner.i_grid[1] + j) % grid_n[1];
        af::tiny<MapFloatType, 4> ppp(0.0);
        for (int k = -1; k < 3; k++) {
          iv_t w = (corner.i_grid[2] + k) % grid_n[2];
          ppp[k+1] = map(u,v,w);
        }
        pp[j+1] = interpolate_at_point(ppp, xn[2]);
      }
      p[i+1] = interpolate_at_point(pp, xn[1]);
    }
    MapFloatType result = interpolate_at_point(p, xn[0]);
    return result;
  }

  template <typename FloatType>
  typename af::c_grid_padded<3>::index_type
  closest_grid_point(
    af::flex_grid<> const& map_accessor,
    fractional<FloatType> const& x_frac)
  {
    af::c_grid_padded<3> c_grid(map_accessor);
    typedef af::c_grid_padded<3>::index_type index_t;
    index_t const& grid_n = c_grid.focus();
    return get_corner<index_t, FloatType>(grid_n, x_frac)
      .closest_grid_point(grid_n);
  }

  template <typename FloatType>
  FloatType
  non_crystallographic_eight_point_interpolation(
    af::const_ref<FloatType, af::flex_grid<> > const& map,
    scitbx::mat3<FloatType> const& gridding_matrix,
    scitbx::vec3<FloatType> const& site_cart,
    bool allow_out_of_bounds=false,
    FloatType const& out_of_bounds_substitute_value=0)
  {
    CCTBX_ASSERT(map.accessor().nd() == 3);
    typedef typename af::flex_grid<>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t map_index(3, 0);
    get_corner<index_t, FloatType> corner(gridding_matrix, site_cart);
    for(unsigned i=0;i<3;i++) {
      if(corner.i_grid[i] < map.accessor().origin()[i] ||
         corner.i_grid[i] >= map.accessor().focus()[i]-1) {
        if(!allow_out_of_bounds) {
          throw error("non_crystallographic_eight_point_interpolation:"
                      " point required for interpolation is out of bounds.");
        }
        else {
          return out_of_bounds_substitute_value;
        }
      }
    }
    FloatType result = 0;
    for(iv_t s0=0;s0<2;s0++) { map_index[0] = (corner.i_grid[0] + s0);
    for(iv_t s1=0;s1<2;s1++) { map_index[1] = (corner.i_grid[1] + s1);
    for(iv_t s2=0;s2<2;s2++) { map_index[2] = (corner.i_grid[2] + s2);
      result += map(map_index) * corner.weight(s0,s1,s2);
    }}}
    return result;
  }

  template <typename FloatType>
  FloatType
  asu_eight_point_interpolation(
    af::const_ref<FloatType, af::flex_grid<> > const& map,
    crystal::direct_space_asu::asu_mappings<FloatType> & am,
    fractional<FloatType> const& site_cart)
  {
    CCTBX_ASSERT(map.accessor().nd() == 3);
    typedef typename af::flex_grid<>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t map_index(3, 0);
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, FloatType> corner(am, grid_n, site_cart);
    FloatType epsilon = scitbx::math::floating_point_epsilon<FloatType>::get() * 10;
    FloatType result = 0;
    for(iv_t s0=0;s0<2;s0++) {
      map_index[0] = (corner.i_grid[0] + s0);
      for(iv_t s1=0;s1<2;s1++) {
        map_index[1] = (corner.i_grid[1] + s1);
        for(iv_t s2=0;s2<2;s2++) {
          map_index[2] = (corner.i_grid[2] + s2);
          // I know comments just *kill* the readers of this file, but I think this needs
          // a bit of explanation:
          // (1) the ASU has grid points which are on the "open face" of the ASU
          // (2) these grid-points are symmetrically related to other, more special grid-points
          // (3) that means that for some grid-points we have to see if they are IN the ASU
          //     then if they are NOT, we have to map them back into the ASU
          // (4) note that there are some grid-points which are "valid" accesses are not
          //     ACTUALLY valid, that's why we call is_valid_index and cross our flindas
          // -j & E. June 22, 2005
          if ( ! map.accessor().is_valid_index(map_index) ) {
            fractional<FloatType> lmap;
            for ( std::size_t i=0; i<3; ++i ) {
              lmap[i] = static_cast<FloatType>(map_index[i]) / grid_n[i];
            }
            cartesian<FloatType> const & xmap = am.process(lmap).mappings().back()[0].mapped_site();
            fractional<FloatType> nxmap = am.unit_cell().fractionalize(xmap);
            for ( std::size_t i=0; i<3; ++i ) {
              if ( std::abs(nxmap[i]) < epsilon )
                nxmap[i] = 0;
              map_index[i] = scitbx::math::float_int_conversions<FloatType, long>::
                                ifloor(nxmap[i] * static_cast<FloatType>(grid_n[i]));
            }
          }
          result += map(map_index) * corner.weight(s0,s1,s2);
        }
      }
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_EIGHT_POINT_INTERPOLATION_H
