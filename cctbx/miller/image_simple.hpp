#ifndef CCTBX_MILLER_IMAGE_SIMPLE_H
#define CCTBX_MILLER_IMAGE_SIMPLE_H

#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <scitbx/mat3.h>
#include <scitbx/vec2.h>
#include <scitbx/array_family/versa.h>
#include <tbxx/error_utils.hpp>

namespace cctbx { namespace miller {

  struct image_simple
  {
    bool store_spots;
    bool set_pixels;
    af::shared<scitbx::vec3<double> > spots;
    af::versa<int, af::flex_grid<> > pixels;

    image_simple(
      bool store_spots_,
      bool set_pixels_)
    :
      store_spots(store_spots_),
      set_pixels(set_pixels_)
    {
      TBXX_ASSERT(store_spots || set_pixels);
    }

    image_simple&
    compute(
      uctbx::unit_cell const& unit_cell,
      af::const_ref<index<> > const& miller_indices,
      af::const_ref<double> const& spot_intensity_factors,
      scitbx::mat3<double> const& crystal_rotation_matrix,
      double ewald_radius,
      double ewald_proximity,
      int signal_max,
      double detector_distance,
      scitbx::vec2<double> detector_size,
      scitbx::vec2<int> detector_pixels,
      unsigned point_spread,
      double gaussian_falloff_scale)
    {
      if (spot_intensity_factors.size() != 0) {
        TBXX_ASSERT(spot_intensity_factors.size() == miller_indices.size());
        TBXX_ASSERT(spot_intensity_factors.all_ge(0));
        TBXX_ASSERT(spot_intensity_factors.all_le(1));
      }
      TBXX_ASSERT(ewald_radius > 0);
      TBXX_ASSERT(detector_size.const_ref().all_gt(0));
      TBXX_ASSERT(detector_pixels.const_ref().all_gt(0));
      TBXX_ASSERT(point_spread > 0);
      TBXX_ASSERT(gaussian_falloff_scale >= 0);
      int dpx = detector_pixels[0];
      int dpy = detector_pixels[1];
      if (set_pixels) {
        pixels.resize(af::flex_grid<>(dpx, dpy), 0);
      }
      int* pixels_beg = pixels.begin();
      double dsx = detector_size[0];
      double dsy = detector_size[1];
      unsigned point_spread_half = point_spread / 2;
      bool point_spread_is_even_value = (point_spread % 2 == 0);
      double circle_radius_sq = point_spread * std::max(dsx/dpx, dsy/dpy) / 2;
      circle_radius_sq *= circle_radius_sq;
      TBXX_ASSERT(circle_radius_sq != 0);
      bool apply_proximity_factor = (ewald_proximity > 0);
      if (!apply_proximity_factor) ewald_proximity *= -1;
      typedef scitbx::vec3<double> v3d;
      for(std::size_t ih=0;ih<miller_indices.size();ih++) {
        v3d rv = unit_cell.reciprocal_space_vector(miller_indices[ih]);
        v3d rvre = crystal_rotation_matrix * rv;
        rvre[2] += ewald_radius; // direct beam anti-parallel (0,0,1)
        double rvre_len = rvre.length();
        double rvre_proximity = std::abs(1 - rvre_len / ewald_radius);
        if (rvre_proximity < ewald_proximity) {
          // http://en.wikipedia.org/wiki/Line-plane_intersection
          if (rvre[2] > 0) {
            double d = -detector_distance / rvre[2];
            double dx = rvre[0] * d;
            double dy = rvre[1] * d;
            if (   std::abs(dx) <= dsx/2
                && std::abs(dy) <= dsy/2) {
              using scitbx::math::ifloor;
              double pxf = (dx/dsx + 0.5) * dpx;
              double pyf = (dy/dsy + 0.5) * dpy;
              if (store_spots) {
                spots.push_back(scitbx::vec3<double>(pxf, pyf, 0));
              }
              if (!set_pixels) continue;
              int pxi = ifloor(pxf);
              int pyi = ifloor(pyf);
              int pxb = pxi - point_spread_half;
              int pyb = pyi - point_spread_half;
              if (point_spread_is_even_value) {
                if (pxf - pxi > 0.5) pxb++;
                if (pyf - pyi > 0.5) pyb++;
              }
              double proximity_factor = 1;
              if (apply_proximity_factor) {
                proximity_factor -= scitbx::fn::pow2(
                  rvre_proximity / ewald_proximity);
                if (proximity_factor <= 0) continue;
              }
              double signal_at_center =
                  signal_max
                * (spot_intensity_factors.size() == 0 ? 1 :
                   spot_intensity_factors[ih])
                * proximity_factor;
              int signal = static_cast<int>(signal_at_center + 0.5);
              double gauss_arg_term = -gaussian_falloff_scale
                                    / circle_radius_sq;
              for(int i=0;i<=point_spread;i++) {
                int pi = pxb + i;
                if (pi < 0 || pi >= dpx) continue;
                int pi0 = pi * dpy;
                for(int j=0;j<=point_spread;j++) {
                  int pj = pyb + j;
                  if (pj < 0 || pj >= dpy) continue;
                  if (point_spread > 2) {
                    double pcx = ((pi + 0.5) / dpx - 0.5) * dsx - dx;
                    double pcy = ((pj + 0.5) / dpy - 0.5) * dsy - dy;
                    double pc_sq = pcx*pcx + pcy*pcy;
                    if (pc_sq > circle_radius_sq) continue;
                    if (gaussian_falloff_scale != 0) {
                      double falloff_factor = std::exp(pc_sq * gauss_arg_term);
                      signal = static_cast<int>(
                        signal_at_center * falloff_factor + 0.5);
                    }
                  }
                  pixels[pi0+pj] = signal;
                }
              }
            }
          }
        }
      }
      return *this;
    }
  };

}} // namespace cctbx::miller

#endif // GUARD
