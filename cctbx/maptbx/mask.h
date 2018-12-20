#ifndef CCTBX_MAPTBX_MASK_H
#define CCTBX_MAPTBX_MASK_H

#include <cstddef>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/uctbx.h>

namespace cctbx { namespace maptbx {

/*
Less efficient (than asu mask) implementation of mask calculation
*/
namespace sm=scitbx::math;
template <typename FloatType>
af::versa<FloatType, af::c_grid<3> > mask(
  af::const_ref<scitbx::vec3<FloatType> > const& sites_frac,
  uctbx::unit_cell const& unit_cell,
  af::tiny<int, 3> const& n_real,
  FloatType const& mask_value_inside_molecule,
  FloatType const& mask_value_outside_molecule,
  af::const_ref<FloatType> const& radii,
  bool const& wrapping
      )
{
  int nx = n_real[0];
  int ny = n_real[1];
  int nz = n_real[2];
  af::versa<FloatType, af::c_grid<3> >
    result(af::c_grid<3>(nx,ny,nz), mask_value_outside_molecule);
  af::ref<FloatType, af::c_grid<3> > result_ref = result.ref();
  af::tiny<FloatType, 6> ucp = unit_cell.parameters();
  FloatType ucs = unit_cell.volume() / (ucp[0]*ucp[1]*ucp[2]);
  for(int j = 0; j < sites_frac.size(); j++) {
    cctbx::fractional<> site_frac = sites_frac[j];
    af::tiny<int, 3> box_min, box_max;
    for(int i = 0; i <= 2; i++) {
      FloatType rf=radii[j]/ucp[i]/(ucs/std::sin(scitbx::deg_as_rad(ucp[i+3])));
      if (wrapping) {
        box_min[i] = sm::nearest_integer(n_real[i]*(site_frac[i]-rf));
        box_max[i] = sm::nearest_integer(n_real[i]*(site_frac[i]+rf));
      } else {
        box_min[i] = std::max(0,
           sm::nearest_integer(n_real[i]*(site_frac[i]-rf)));
        box_max[i] = std::min(n_real[i]-1,
           sm::nearest_integer(n_real[i]*(site_frac[i]+rf)));
      }

    }
    for(int kx = box_min[0]; kx < box_max[0]; kx++) {
      FloatType xn=FloatType(kx)/n_real[0];
      int kx_ = sm::mod_positive(kx, n_real[0]);
      for(int ky = box_min[1]; ky < box_max[1]; ky++) {
        FloatType yn=FloatType(ky)/n_real[1];
        int ky_ = sm::mod_positive(ky, n_real[1]);
        for(int kz = box_min[2]; kz < box_max[2]; kz++) {
          FloatType zn=FloatType(kz)/n_real[2];
          FloatType distance = unit_cell.distance(site_frac,
            cctbx::fractional<>(xn,yn,zn));
          int kz_ = sm::mod_positive(kz, n_real[2]);
          if(distance < radii[j]) {
            result_ref(kx_,ky_,kz_)=mask_value_inside_molecule;
          }
    }}}
  }
  return result;
}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_UTILS_H
