#ifndef CCTBX_MAPTBX_INTERPOLATION_H
#define CCTBX_MAPTBX_INTERPOLATION_H

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
      FloatType weights_[3][2];
  };

  template <
    typename MapFloatType,
    typename SiteFloatType>
  af::tiny<MapFloatType, 4>
  eight_point_interpolation_with_gradients(
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac,
    scitbx::vec3<SiteFloatType> const& step)
  {
    typedef af::c_grid_padded<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    MapFloatType result = 0;
    MapFloatType f_000, f_100, f_010, f_110, f_001, f_101, f_011, f_111;
    for(iv_t s0=0;s0<2;s0++) { iv_t i0 = (corner.i_grid[0] + s0) % grid_n[0];
    for(iv_t s1=0;s1<2;s1++) { iv_t i1 = (corner.i_grid[1] + s1) % grid_n[1];
    for(iv_t s2=0;s2<2;s2++) { iv_t i2 = (corner.i_grid[2] + s2) % grid_n[2];
      MapFloatType map_value = map(i0,i1,i2);
      result += map_value * corner.weight(s0,s1,s2);
      if(s0==0&&s1==0&&s2==0) f_000 = map_value;
      if(s0==1&&s1==0&&s2==0) f_100 = map_value;
      if(s0==0&&s1==1&&s2==0) f_010 = map_value;
      if(s0==1&&s1==1&&s2==0) f_110 = map_value;
      if(s0==0&&s1==0&&s2==1) f_001 = map_value;
      if(s0==1&&s1==0&&s2==1) f_101 = map_value;
      if(s0==0&&s1==1&&s2==1) f_011 = map_value;
      if(s0==1&&s1==1&&s2==1) f_111 = map_value;
    }}}
    MapFloatType x = corner.weights_[0][1];
    MapFloatType y = corner.weights_[1][1];
    MapFloatType z = corner.weights_[2][1];
    MapFloatType f_x00 = (1-x)*f_000 + x*f_100;
    MapFloatType f_x01 = (1-x)*f_001 + x*f_101;
    MapFloatType f_0y0 = (1-y)*f_000 + y*f_010;
    MapFloatType f_0y1 = (1-y)*f_001 + y*f_011;
    MapFloatType f_x10 = (1-x)*f_010 + x*f_110;
    MapFloatType f_x11 = (1-x)*f_011 + x*f_111;
    MapFloatType f_1y0 = (1-y)*f_100 + y*f_110;
    MapFloatType f_1y1 = (1-y)*f_101 + y*f_111;
    MapFloatType f_xy0 = (1-y)*f_x00 + y*f_x10;
    MapFloatType f_x0z = (1-z)*f_x00 + z*f_x01;
    MapFloatType f_0yz = (1-z)*f_0y0 + z*f_0y1;
    MapFloatType f_xy1 = (1-y)*f_x01 + y*f_x11;
    MapFloatType f_x1z = (1-z)*f_x10 + z*f_x11;
    MapFloatType f_1yz = (1-z)*f_1y0 + z*f_1y1;
    // Comment out for performance
    CCTBX_ASSERT( std::abs((1-z)*f_xy0+z*f_xy1-result)<1.e-6 );
    CCTBX_ASSERT( std::abs((1-x)*f_0yz+x*f_1yz-result)<1.e-6 );
    CCTBX_ASSERT( std::abs((1-y)*f_x0z+y*f_x1z-result)<1.e-6 );
    MapFloatType gx = (f_1yz-f_0yz) / step[0];
    MapFloatType gy = (f_x1z-f_x0z) / step[1];
    MapFloatType gz = (f_xy1-f_xy0) / step[2];
    return af::tiny<MapFloatType, 4>(result, gx,gy,gz);
  }

  template <
    typename MapFloatType,
    typename SiteFloatType>
  af::tiny<MapFloatType, 4>
  quadratic_interpolation_with_gradients(
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& map,
    //af::const_ref<MapFloatType, af::flex_grid<> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac,
    scitbx::vec3<SiteFloatType> const& step)
  {
    typedef af::c_grid_padded<3>::index_type index_t;
    //typedef af::flex_grid<>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    MapFloatType f_000, f_100, f_010, f_110, f_001, f_101, f_011, f_111;
    MapFloatType f_m100, f_0m10, f_00m1;
    MapFloatType f_200, f_020, f_002;
    for(int s0=-1;s0<3;s0++) { int i0 = (corner.i_grid[0] + s0) % grid_n[0];
    for(int s1=-1;s1<3;s1++) { int i1 = (corner.i_grid[1] + s1) % grid_n[1];
    for(int s2=-1;s2<3;s2++) { int i2 = (corner.i_grid[2] + s2) % grid_n[2];
      MapFloatType map_value = map(i0,i1,i2);
      if(s0==0&&s1==0&&s2==0) f_000 = map_value;
      if(s0==1&&s1==0&&s2==0) f_100 = map_value;
      if(s0==0&&s1==1&&s2==0) f_010 = map_value;
      if(s0==1&&s1==1&&s2==0) f_110 = map_value;
      if(s0==0&&s1==0&&s2==1) f_001 = map_value;
      if(s0==1&&s1==0&&s2==1) f_101 = map_value;
      if(s0==0&&s1==1&&s2==1) f_011 = map_value;
      if(s0==1&&s1==1&&s2==1) f_111 = map_value;
      if(s0==-1&&s1== 0&&s2== 0) f_m100 = map_value;
      if(s0== 0&&s1==-1&&s2== 0) f_0m10 = map_value;
      if(s0== 0&&s1== 0&&s2==-1) f_00m1 = map_value;
      if(s0==2&&s1==0&&s2==0) f_200 = map_value;
      if(s0==0&&s1==2&&s2==0) f_020 = map_value;
      if(s0==0&&s1==0&&s2==2) f_002 = map_value;
    }}}
    MapFloatType x = corner.weights_[0][1];
    MapFloatType y = corner.weights_[1][1];
    MapFloatType z = corner.weights_[2][1];
    MapFloatType t = (f_111+f_100+f_010+f_001-f_011-f_110-f_101-f_000) / 8;
    MapFloatType f_000_bar = f_000-t;
    MapFloatType f_100_bar = f_100-t;
    MapFloatType f_010_bar = f_010-t;
    MapFloatType f_110_bar = f_110-t;
    MapFloatType f_001_bar = f_001-t;
    MapFloatType f_101_bar = f_101-t;
    MapFloatType f_011_bar = f_011-t;
    MapFloatType f_111_bar = f_111-t;
    MapFloatType axy = (f_110_bar+f_000_bar-f_100_bar-f_010_bar)/2;
    MapFloatType ayz = (f_011_bar+f_000_bar-f_010_bar-f_001_bar)/2;
    MapFloatType axz = (f_101_bar+f_000_bar-f_100_bar-f_001_bar)/2;
    MapFloatType c = f_000_bar;
    MapFloatType axx,ayy,azz, bx,by,bz;
    if(x>=0 && x<0.5) {
      axx = (f_100_bar-2*f_000_bar+f_m100)/2;
      bx  = (f_100_bar-f_m100)/2;
    }
    else {
      axx = (f_200-2*f_100_bar+f_000_bar)/2;
      bx  = (4*f_100_bar-f_200-3*f_000_bar)/2;
    }
    if(y>=0 && y<0.5) {
      ayy = (f_010_bar-2*f_000+f_0m10)/2;
      by  = (f_010_bar-f_0m10)/2;
    }
    else {
      ayy = (f_020-2*f_010_bar+f_000_bar)/2;
      by  = (4*f_010_bar-f_020-3*f_000_bar)/2;
    }
    if(z>=0 && z<0.5) {
      azz = (f_001_bar-2*f_000_bar+f_00m1)/2.;
      bz  = (f_001_bar-f_00m1)/2;
    }
    else {
      azz = (f_002-2*f_001_bar+f_000_bar)/2;
      bz  = (4*f_001_bar-f_002-3*f_000_bar)/2;
    }
    MapFloatType result = axx*x*x +
                          ayy*y*y +
                          azz*z*z +
                          2*axy*x*y +
                          2*axz*x*z +
                          2*ayz*y*z +
                          bx*x + by*y + bz*z + c;
    MapFloatType gx = (2*x*axx + 2*y*axy + 2*z*axz + bx) / step[0];
    MapFloatType gy = (2*x*axy + 2*y*ayy + 2*z*ayz + by) / step[1];
    MapFloatType gz = (2*x*axz + 2*y*ayz + 2*z*azz + bz) / step[2];
    return af::tiny<MapFloatType, 4>(result, gx,gy,gz);
  }

  template <
    typename MapFloatType,
    typename SiteFloatType>
  MapFloatType
  eight_point_interpolation(
    af::const_ref<MapFloatType, af::flex_grid<> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac)
  {
    typedef af::flex_grid<>::index_type index_t;
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

  template <typename FloatType>
  FloatType
  cubic(FloatType t, FloatType fm1, FloatType f0, FloatType f1, FloatType f2) {
    FloatType a0 = f0;
    FloatType a1 = (-1*f2+6*f1-3*f0-2*fm1)/6.;
    FloatType a2 = (f1-2*f0+fm1)/2.;
    FloatType a3 = (f2-3*f1+3*f0-fm1)/6.;
    //
    // Other option
    //FloatType a0 = f0;
    //FloatType a1 = (f1-fm1)/2.;
    //FloatType a2 = (-1.*f2+4*f1-5*f0+2*fm1)/2.;
    //FloatType a3 = (f2-3*f1+3*f0-fm1)/2.;
    return a0 + t*(a1 + t*(a2 + a3*t));
  }

  template <typename FloatType>
  FloatType
  gcubic(FloatType t, FloatType fm1, FloatType f0, FloatType f1, FloatType f2) {
    FloatType a1 = (-1*f2+6*f1-3*f0-2*fm1)/6.;
    FloatType a2 = (f1-2*f0+fm1)/2.;
    FloatType a3 = (f2-3*f1+3*f0-fm1)/6.;
    //
    // Other option
    //FloatType a1 = (f1-fm1)/2.;
    //FloatType a2 = (-1.*f2+4*f1-5*f0+2*fm1)/2.;
    //FloatType a3 = (f2-3*f1+3*f0-fm1)/2.;
    return a1 + t*(2*a2 + 3*a3*t);
  }

  template <typename FloatType>
  FloatType
  fxpq(FloatType f[4][4][4], FloatType x, int p, int q) {
    p=p+1;
    q=q+1;
    return cubic(x, f[-1+1][p][q], f[0+1][p][q], f[1+1][p][q], f[2+1][p][q]);
  }

  template <typename FloatType>
  FloatType
  fqyp(FloatType f[4][4][4], FloatType y, int p, int q) {
    p=p+1;
    q=q+1;
    return cubic(y, f[q][-1+1][p], f[q][0+1][p], f[q][1+1][p], f[q][2+1][p]);
  }

  template <typename FloatType>
  FloatType
  fpqz(FloatType f[4][4][4], FloatType z, int p, int q) {
    p=p+1;
    q=q+1;
    return cubic(z, f[p][q][-1+1], f[p][q][0+1], f[p][q][1+1], f[p][q][2+1]);
  }

  template <typename FloatType>
  FloatType
  fxyq(FloatType f[4][4][4], FloatType x, FloatType y, int q) {
    return cubic(
      y,
      fxpq(f, x, -1, q),
      fxpq(f, x,  0, q),
      fxpq(f, x,  1, q),
      fxpq(f, x,  2, q));
  }

  template <typename FloatType>
  FloatType
  fqyz(FloatType f[4][4][4], FloatType y, FloatType z, int q) {
    return cubic(
      z,
      fqyp(f, y, -1, q),
      fqyp(f, y,  0, q),
      fqyp(f, y,  1, q),
      fqyp(f, y,  2, q));
  }

  template <typename FloatType>
  FloatType
  fxqz(FloatType f[4][4][4], FloatType x, FloatType z, int q) {
    return cubic(
      x,
      fpqz(f, z, -1, q),
      fpqz(f, z,  0, q),
      fpqz(f, z,  1, q),
      fpqz(f, z,  2, q));
  }

  template <
  typename MapFloatType,
  typename SiteFloatType>
  af::tiny<MapFloatType, 4>
  tricubic_interpolation_with_gradients(
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac,
    scitbx::vec3<SiteFloatType> const& step)
  {
    using namespace std;
    typedef af::c_grid_padded<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor().focus();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    MapFloatType f[4][4][4];
    for (int i = -1; i < 3; i++) {
      iv_t u = (corner.i_grid[0] + i) % grid_n[0];
      for (int j = -1; j < 3; j++) {
        iv_t v = (corner.i_grid[1] + j) % grid_n[1];
        for (int k = -1; k < 3; k++) {
          iv_t w = (corner.i_grid[2] + k) % grid_n[2];
          f[i+1][j+1][k+1] = map(u,v,w);
        }
      }
    }
    SiteFloatType x = corner.weights_[0][1];
    SiteFloatType y = corner.weights_[1][1];
    SiteFloatType z = corner.weights_[2][1];
    // All three must be the same
    // Comment out any two for performance
    MapFloatType r0 = cubic<MapFloatType>(
      z, fxyq(f, x,y,-1), fxyq(f, x,y,0),
         fxyq(f, x,y,1),  fxyq(f, x,y,2));
    //MapFloatType r1 = cubic<MapFloatType>(
    //  x, fqyz(f, y,z,-1), fqyz(f, y,z,0),
    //     fqyz(f, y,z,1),  fqyz(f, y,z,2));
    //MapFloatType r2 = cubic<MapFloatType>(
    //  y, fxqz(f, x,z,-1), fxqz(f, x,z,0),
    //     fxqz(f, x,z,1),  fxqz(f, x,z,2));
    //CCTBX_ASSERT(std::abs(r0-r2)<1.e-6);
    //CCTBX_ASSERT(std::abs(r0-r1)<1.e-6);
    MapFloatType gx = gcubic<MapFloatType>(
      x, fqyz(f, y,z,-1), fqyz(f, y,z,0),
         fqyz(f, y,z,1),  fqyz(f, y,z,2));
    MapFloatType gy = gcubic<MapFloatType>(
      y, fxqz(f, x,z,-1), fxqz(f, x,z,0),
         fxqz(f, x,z,1),  fxqz(f, x,z,2));
    MapFloatType gz = gcubic<MapFloatType>(
      z, fxyq(f, x,y,-1), fxyq(f, x,y,0),
         fxyq(f, x,y,1),  fxyq(f, x,y,2));
    return af::tiny<MapFloatType, 4>(r0, gx/step[0],gy/step[1],gz/step[2]);
  }

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
        xn[k] = fmod((1.+x_frac[k])*static_cast<SiteFloatType>(grid_n[k]), 1.0);
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

  template <
    typename MapFloatType,
    typename SiteFloatType>
  MapFloatType
  tricubic_interpolation(
    af::const_ref<MapFloatType, af::c_grid<3> > const& map,
    scitbx::vec3<SiteFloatType> const& x_frac)
  {
    using scitbx::math::interpolate_at_point;
    using namespace std;
    typedef af::c_grid<3>::index_type index_t;
    typedef typename index_t::value_type iv_t;
    index_t const& grid_n = map.accessor();
    get_corner<index_t, SiteFloatType> corner(grid_n, x_frac);
    af::tiny<MapFloatType, 4> p(0.0);
    af::tiny<SiteFloatType, 3> xn;
    for (unsigned k = 0; k < 3; k++) {
      if (x_frac[k] < 0) {
        xn[k] = fmod((1.+x_frac[k])*static_cast<SiteFloatType>(grid_n[k]), 1.0);
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

#endif // CCTBX_MAPTBX_INTERPOLATION_H
