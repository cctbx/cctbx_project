#ifndef CCTBX_MAPTBX_REAL_SPACE_GRADIENTS_SIMPLE_2_H
#define CCTBX_MAPTBX_REAL_SPACE_GRADIENTS_SIMPLE_2_H

#include <cctbx/maptbx/eight_point_interpolation.h>

namespace cctbx { namespace maptbx {

  template <
    typename FloatType >
  FloatType
  real_space_target_simple_2(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map_target,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& map_current,
    FloatType const& box_size,
    af::const_ref<scitbx::vec3<FloatType> > const& sites_frac)
  {
    int nx = static_cast<int>(map_target.accessor().focus()[0]);
    int ny = static_cast<int>(map_target.accessor().focus()[1]);
    int nz = static_cast<int>(map_target.accessor().focus()[2]);
    FloatType rp[3];
    for(unsigned i=0;i<3;i++) {
      rp[i] = static_cast<FloatType>(unit_cell.reciprocal_parameters()[i]);
    }
    FloatType result = 0;
    for(std::size_t i_site=0;i_site<sites_frac.size();i_site++) {
      FloatType coas = box_size*rp[0];
      FloatType cobs = box_size*rp[1];
      FloatType cocs = box_size*rp[2];
      cctbx::fractional<> const& site_frac = sites_frac[i_site];
      FloatType xfi= static_cast<FloatType>(site_frac[0]);
      FloatType yfi= static_cast<FloatType>(site_frac[1]);
      FloatType zfi= static_cast<FloatType>(site_frac[2]);
      //int x1box=0;
      //int x2box=nx-1;
      //int y1box=0;
      //int y2box=ny-1;
      //int z1box=0;
      //int z2box=nz-1;

      std::cout<<""<<std::endl;
      std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
      int x1box=std::max(0,   scitbx::math::float_int_conversions<FloatType, int>::ifloor(nx*(xfi-coas)));
      int x2box=std::min(nx-1,scitbx::math::float_int_conversions<FloatType, int>::iceil( nx*(xfi+coas)));
      int y1box=std::max(0,   scitbx::math::float_int_conversions<FloatType, int>::ifloor(ny*(yfi-cobs)));
      int y2box=std::min(ny-1,scitbx::math::float_int_conversions<FloatType, int>::iceil( ny*(yfi+cobs)));
      int z1box=std::max(0,   scitbx::math::float_int_conversions<FloatType, int>::ifloor(nz*(zfi-cocs)));
      int z2box=std::min(nz-1,scitbx::math::float_int_conversions<FloatType, int>::iceil( nz*(zfi+cocs)));
      std::cout<<x1box<<" "<<y1box<<" "<<z1box<<std::endl;
      std::cout<<x2box<<" "<<y2box<<" "<<z2box<<std::endl;
      for(int kx = x1box; kx <= x2box; kx++) {
        int mx = scitbx::math::mod_positive(kx,nx);
        int mxny = mx*ny;
        for(int ky = y1box; ky <= y2box; ky++) {
          int my = scitbx::math::mod_positive(ky,ny);
          int mxnypmynz = (mxny+my)*nz;
          for (int kz = z1box; kz <= z2box; kz++) {
            int mz = scitbx::math::mod_positive(kz,nz);
            FloatType diff = map_target[mxnypmynz+mz]-map_current[mxnypmynz+mz];
            result += diff * diff; // (mx*NY+my)*NZ+mz
          }
        }
      }
    }
    return result/sites_frac.size();
  }

  template <
    typename MapFloatType,
    typename SiteFloatType>
  af::shared<scitbx::vec3<SiteFloatType> >
  real_space_gradients_simple_2(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
    af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
    SiteFloatType delta)
  {
    CCTBX_ASSERT(delta > 0);
    typedef scitbx::vec3<SiteFloatType> v3_t;
    af::shared<v3_t> result(sites_cart.size(), af::init_functor_null<v3_t>());
    v3_t* res = result.begin();
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++,res++) {
      v3_t piv = sites_cart[i_site];
      v3_t piv_d = piv;
      for(unsigned i_axis=0;i_axis<3;i_axis++) {
        MapFloatType densities[2];
        for(unsigned i_sign=0;i_sign<2;i_sign++) {
          piv_d[i_axis] = (i_sign == 0 ? piv[i_axis] + delta
                                       : piv[i_axis] - delta);
          fractional<SiteFloatType> site_frac = unit_cell.fractionalize(piv_d);
          densities[i_sign] = eight_point_interpolation(
            density_map, site_frac);
        }
        piv_d[i_axis] = piv[i_axis];
        (*res)[i_axis] = (densities[0] - densities[1]) / (2 * delta);
      }
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // GUARD
