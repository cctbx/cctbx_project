#ifndef CCTBX_MAPTBX_REAL_SPACE_GRADIENTS_SIMPLE_H
#define CCTBX_MAPTBX_REAL_SPACE_GRADIENTS_SIMPLE_H

#include <cctbx/maptbx/eight_point_interpolation.h>

namespace cctbx { namespace maptbx {

  template <
    typename MapFloatType,
    typename SiteFloatType>
  MapFloatType
  real_space_target_simple(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
    af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
    af::const_ref<bool> const& selection)
  {
    MapFloatType result = 0;
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      if(selection[i_site]) {
        result += eight_point_interpolation(
          density_map,
          unit_cell.fractionalize(sites_cart[i_site]));
      }
    }
    return result;
  }

  template <
    typename MapFloatType,
    typename SiteFloatType>
  af::shared<MapFloatType>
  real_space_target_simple_per_site(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
    af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart)
  {
    af::shared<MapFloatType> result(
      sites_cart.size(), af::init_functor_null<MapFloatType>());
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++) {
      result[i_site] = eight_point_interpolation(
        density_map,
        unit_cell.fractionalize(sites_cart[i_site]));
    }
    return result;
  }

  template <
    typename MapFloatType,
    typename SiteFloatType>
  af::shared<scitbx::vec3<SiteFloatType> >
  real_space_gradients_simple(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<MapFloatType, af::c_grid_padded<3> > const& density_map,
    af::const_ref<scitbx::vec3<SiteFloatType> > const& sites_cart,
    SiteFloatType delta,
    af::const_ref<bool> const& selection)
  {
    CCTBX_ASSERT(delta > 0);
    typedef scitbx::vec3<SiteFloatType> v3_t;
    af::shared<v3_t> result(sites_cart.size(), af::init_functor_null<v3_t>());
    v3_t* res = result.begin();
    for(std::size_t i_site=0;i_site<sites_cart.size();i_site++,res++) {
      if(selection[i_site]) {
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
    }
    return result;
  }

}} // namespace cctbx::maptbx

#endif // GUARD
