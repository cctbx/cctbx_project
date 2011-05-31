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
    scitbx::vec2<int> detector_pixels)
  {
    TBXX_ASSERT(ewald_radius > 0);
    TBXX_ASSERT(detector_pixels.const_ref().all_gt(0));
    int dpx = detector_pixels[0];
    int dpy = detector_pixels[1];
    std::string result(dpx * dpy * 3, static_cast<char>(255U));
    double dsx = detector_size[0];
    double dsy = detector_size[1];
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
            using scitbx::math::iround;
            int px = std::min(iround((dx/dsx + 1/2.) * dpx), dpx-1);
            int py = std::min(iround((dy/dsy + 1/2.) * dpy), dpy-1);
            for(int i=-1;i<=1;i++) {
              int pi = px + i;
              if (pi < 0 || pi >= dpx) continue;
              for(int j=-1;j<=1;j++) {
                int pj = py + j;
                if (pj < 0 || pj >= dpy) continue;
                int pijk = (pi*dpy+pj)*3;
                for(int k=0;k<3;k++,pijk++) {
                  result[pijk] = static_cast<char>(0U);
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
