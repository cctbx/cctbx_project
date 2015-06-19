#ifndef CCTBX_MAPTBX_SPHERICITY_H
#define CCTBX_MAPTBX_SPHERICITY_H

#include <cctbx/uctbx.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/matrix/eigensystem.h>

namespace cctbx { namespace maptbx {

namespace sm=scitbx::math;
namespace es=scitbx::matrix::eigensystem;

//! Compute peak sphericity as the ratio of min/max eigen-values of spericity
//! tensor using coordinates of grid points in a box around peak center selected
//! above one-nth (n=3) of peak max height value. Tensor elements are weighted
//! with map values at respective grid nodes though it isn't clear whether it is
//! necessary conceptially or useful practically.
//! Known (fixable) limitation: may consider two closely standing peaks as one.

template <typename FloatType>
scitbx::sym_mat3<FloatType>
  sphericity_tensor(
    af::const_ref<FloatType, af::c_grid<3> > const& map_data,
    uctbx::unit_cell const& unit_cell,
    FloatType const& radius,
    cctbx::fractional<> const& site_frac)
{
  scitbx::sym_mat3<FloatType> result;
  result.fill(0);
  af::tiny<int, 3> const& n_real = map_data.accessor();
  af::tiny<int, 3> box_min, box_max, n_center;
  af::tiny<FloatType, 6> ucp = unit_cell.parameters();
    FloatType ucs = unit_cell.volume() / (ucp[0]*ucp[1]*ucp[2]);
  for(int i = 0; i <= 2; i++) {
    FloatType rf=radius/ucp[i]/(ucs/std::sin(scitbx::deg_as_rad(ucp[i+3])));
    box_min[i] = sm::nearest_integer(n_real[i]*(site_frac[i]-rf));
    box_max[i] = sm::nearest_integer(n_real[i]*(site_frac[i]+rf));
    n_center[i] = sm::mod_positive(
      sm::nearest_integer(n_real[i]*site_frac[i]), n_real[i]);
  }
  //FloatType min_peak_value = map_data(n_center)/3.;
  for(int kx = box_min[0]; kx < box_max[0]; kx++) {
    FloatType xn=site_frac[0]-FloatType(kx)/n_real[0];
    int mx = sm::mod_positive(kx, n_real[0]);
    for(int ky = box_min[1]; ky < box_max[1]; ky++) {
      FloatType yn=site_frac[1]-FloatType(ky)/n_real[1];
      int my = sm::mod_positive(ky, n_real[1]);
      for(int kz = box_min[2]; kz < box_max[2]; kz++) {
        FloatType zn=site_frac[2]-FloatType(kz)/n_real[2];
        int mz = sm::mod_positive(kz, n_real[2]);
        FloatType map_value = map_data(mx,my,mz);
        cctbx::cartesian<> sc = unit_cell.orthogonalize(
          cctbx::fractional<>(xn,yn,zn));
        //if(map_value>min_peak_value) {
          FloatType w = map_value;
          result(0,0) += w*sc[0]*sc[0];
          result(1,1) += w*sc[1]*sc[1];
          result(2,2) += w*sc[2]*sc[2];
          result(0,1) += w*sc[0]*sc[1];
          result(0,2) += w*sc[0]*sc[2];
          result(1,2) += w*sc[1]*sc[2];
        //}
  }}}
  return result;
}

template <typename FloatType>
af::shared<FloatType>
  sphericity(
    af::const_ref<FloatType, af::c_grid<3> > const& map_data,
    uctbx::unit_cell const& unit_cell,
    FloatType const& radius,
    af::const_ref<scitbx::vec3<FloatType> > const& sites_frac)
{
    af::tiny<int, 3> const& n_real = map_data.accessor();
    af::shared<FloatType> result;
    result.resize(sites_frac.size(), 0.0);
    scitbx::sym_mat3<FloatType> st;
    for(int j = 0; j < sites_frac.size(); j++) {
      cctbx::fractional<> site_frac = sites_frac[j];
      st = sphericity_tensor(map_data, unit_cell, radius, site_frac);
      af::shared<FloatType> ev = es::real_symmetric<FloatType>(st).values();
      FloatType den = af::max(ev.ref());
      if(den != 0) result[j] = af::min(ev.ref())/den;
    }
    return result;
}

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_SPHERICITY_H
