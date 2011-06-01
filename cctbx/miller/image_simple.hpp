#ifndef CCTBX_MILLER_IMAGE_SIMPLE_H
#define CCTBX_MILLER_IMAGE_SIMPLE_H

#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <scitbx/mat3.h>
#include <scitbx/vec2.h>
#include <tbxx/error_utils.hpp>

namespace cctbx { namespace miller {

  inline
  std::string
  image_simple(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<index<> > const& miller_indices,
    scitbx::mat3<double> const& crystal_rotation_matrix,
    double ewald_radius,
    double ewald_proximity,
    double detector_distance,
    scitbx::vec2<double> detector_size,
    scitbx::vec2<int> detector_pixels,
    unsigned point_spread)
  {
    TBXX_ASSERT(ewald_radius > 0);
    TBXX_ASSERT(detector_pixels.const_ref().all_gt(0));
    TBXX_ASSERT(point_spread > 0);
    int dpx = detector_pixels[0];
    int dpy = detector_pixels[1];
    std::string result(dpx * dpy * 3, static_cast<char>(255U));
    double dsx = detector_size[0];
    double dsy = detector_size[1];
    unsigned ps2 = point_spread / 2;
    bool point_spread_is_even_value = (point_spread % 2 == 0);
    typedef scitbx::vec3<double> v3d;
    for(std::size_t ih=0;ih<miller_indices.size();ih++) {
      v3d rv = unit_cell.reciprocal_space_vector(miller_indices[ih]);
      v3d rvre = crystal_rotation_matrix * rv;
      rvre[2] += ewald_radius; // direct beam parallel (0,0,1)
      double rvre_len = rvre.length();
      double rvre_proximity = rvre_len / ewald_radius;
      if (std::abs(1-rvre_proximity) <= ewald_proximity) {
        // http://en.wikipedia.org/wiki/Line-plane_intersection
        if (rvre[2] > 0) {
          double d = detector_distance / rvre[2];
          double dx = rvre[0] * d;
          double dy = rvre[1] * d;
          if (   std::abs(dx) <= dsx/2
              && std::abs(dy) <= dsy/2) {
            using scitbx::math::ifloor;
            double pxf = (dx/dsx + 0.5) * dpx;
            double pyf = (dy/dsy + 0.5) * dpy;
            int pxi = ifloor(pxf);
            int pyi = ifloor(pyf);
            int pxb = pxi - ps2;
            int pyb = pyi - ps2;
            if (point_spread_is_even_value) {
              if (pxf - pxi > 0.5) pxb++;
              if (pyf - pyi > 0.5) pyb++;
            }
            for(int i=0;i<point_spread;i++) {
              int pi = pxb + i;
              if (pi < 0 || pi >= dpx) continue;
              for(int j=0;j<point_spread;j++) {
                int pj = pyb + j;
                if (pj < 0 || pj >= dpy) continue;
                int pij = (pi*dpy+pj)*3;
                for(int k=0;k<3;k++) {
                  result[pij+k] = static_cast<char>(0U);
                }
              }
            }
          }
        }
      }
    }
    return result;
  }

}} // namespace cctbx::miller

#endif // GUARD
