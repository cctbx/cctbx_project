#ifndef CCTBX_MAPTBX_REAL_SPACE_TARGET_AND_GRADIENTS_H
#define CCTBX_MAPTBX_REAL_SPACE_TARGET_AND_GRADIENTS_H

#include <cctbx/maptbx/eight_point_interpolation.h>

namespace cctbx { namespace maptbx {

class target_and_gradients {
public:
  target_and_gradients(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<double, af::c_grid_padded<3> > const& map_target,
    af::const_ref<double, af::c_grid_padded<3> > const& map_current,
    double const& step,
    af::const_ref<scitbx::vec3<double> > const& sites_frac)
  {
    int nx = static_cast<int>(map_target.accessor().focus()[0]);
    int ny = static_cast<int>(map_target.accessor().focus()[1]);
    int nz = static_cast<int>(map_target.accessor().focus()[2]);
    af::versa<double, af::c_grid_padded<3> > diff_density_array;
    diff_density_array.resize(af::c_grid_padded<3>(nx,ny,nz), 0);
    int x1box=0   ;
    int x2box=nx-1;
    int y1box=0   ;
    int y2box=ny-1;
    int z1box=0   ;
    int z2box=nz-1;
    target_ = 0;
    double* diff_density_array_ = diff_density_array.begin();
    for(int kx = x1box; kx <= x2box; kx++) {
      int mx = scitbx::math::mod_positive(kx,nx);
      int mxny = mx*ny;
      for(int ky = y1box; ky <= y2box; ky++) {
        int my = scitbx::math::mod_positive(ky,ny);
        int mxnypmynz = (mxny+my)*nz;
        for(int kz = z1box; kz <= z2box; kz++) {
          int mz = scitbx::math::mod_positive(kz,nz);
          double diff = map_target[mxnypmynz+mz]-map_current[mxnypmynz+mz];
          target_ += (diff * diff); // (mx*NY+my)*NZ+mz
          diff_density_array_[mxnypmynz+mz] = -2 * diff;
        }
      }
    }
    // gradients
    af::const_ref<double, af::c_grid_padded<3> > diffd =
      diff_density_array.const_ref();
    cctbx::cartesian<> step_x = cctbx::cartesian<>(step,0,0);
    cctbx::cartesian<> step_y = cctbx::cartesian<>(0,step,0);
    cctbx::cartesian<> step_z = cctbx::cartesian<>(0,0,step);
    double two_step = 2*step;
    gradients_.resize(sites_frac.size(), scitbx::vec3<double>(0,0,0));
    for(std::size_t i_site=0;i_site<sites_frac.size();i_site++) {
      cctbx::fractional<> const& site_frac = sites_frac[i_site];
      cctbx::cartesian<> site_cart = unit_cell.orthogonalize(site_frac);
      cctbx::fractional<> sxp = unit_cell.fractionalize(site_cart+step_x);
      cctbx::fractional<> sxm = unit_cell.fractionalize(site_cart-step_x);
      cctbx::fractional<> syp = unit_cell.fractionalize(site_cart+step_y);
      cctbx::fractional<> sym = unit_cell.fractionalize(site_cart-step_y);
      cctbx::fractional<> szp = unit_cell.fractionalize(site_cart+step_z);
      cctbx::fractional<> szm = unit_cell.fractionalize(site_cart-step_z);
      double gx = (eight_point_interpolation(diffd, sxp) -
                   eight_point_interpolation(diffd, sxm)) / two_step;
      double gy = (eight_point_interpolation(diffd, syp) -
                   eight_point_interpolation(diffd, sym)) / two_step;
      double gz = (eight_point_interpolation(diffd, szp) -
                   eight_point_interpolation(diffd, szm)) / two_step;
      gradients_[i_site] = scitbx::vec3<double>(gx,gy,gz);
    }
  }

  double target() { return target_; }
  af::shared<scitbx::vec3<double> > gradients() { return gradients_; }

protected:
  double target_;
  af::shared<scitbx::vec3<double> > gradients_;
};


}} // namespace cctbx::maptbx

#endif // GUARD
