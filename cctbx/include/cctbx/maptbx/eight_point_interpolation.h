#ifndef CCTBX_MAPTBX_EIGHT_POINT_INTERPOLATION_H
#define CCTBX_MAPTBX_EIGHT_POINT_INTERPOLATION_H

#include <cctbx/coordinates.h>
#include <cctbx/math/mod.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/math/utils.h>

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
          i_grid[i] = math::mod_positive(
            ixn, static_cast<SignedIntType>(grid_n[i]));
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

  template <typename FloatType>
  FloatType
  eight_point_interpolation(
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map,
    fractional<FloatType> const& x_frac)
  {
    typedef af::c_grid_padded<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, FloatType> corner(grid_n, x_frac);
    FloatType result = 0;
    for(iv_t s0=0;s0<2;s0++) { iv_t i0 = (corner.i_grid[0] + s0) % grid_n[0];
    for(iv_t s1=0;s1<2;s1++) { iv_t i1 = (corner.i_grid[1] + s1) % grid_n[1];
    for(iv_t s2=0;s2<2;s2++) { iv_t i2 = (corner.i_grid[2] + s2) % grid_n[2];
      result += map(i0,i1,i2) * corner.weight(s0,s1,s2);
    }}}
    return result;
  }

  template <typename FloatType>
  typename af::c_grid_padded<3>::index_type
  closest_grid_point(
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map,
    fractional<FloatType> const& x_frac)
  {
    typedef af::c_grid_padded<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    return get_corner<index_t, FloatType>(grid_n, x_frac)
      .closest_grid_point(grid_n);
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_EIGHT_POINT_INTERPOLATION_H
