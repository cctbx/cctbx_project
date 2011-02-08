#ifndef SCITBX_MATH_INTERPOLATION_H
#define SCITBX_MATH_INTERPOLATION_H

#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>

namespace scitbx { namespace math {

  namespace af = scitbx::af;

  // http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  // http://www.mvps.org/directx/articles/catmull
  // XXX: convert to template with arbitary vector dimensionality
  template <typename PointType>
  af::shared< PointType > interpolate_catmull_rom_spline (
    PointType const& p0,
    PointType const& p1,
    PointType const& p2,
    PointType const& p3,
    unsigned n_points)
  {
    SCITBX_ASSERT(n_points >= 1);
    af::shared< PointType > spline(n_points);
    for (unsigned i = 1; i <= n_points; i++) {
      double t = ((double) i) / n_points;
      double t2 = t * t;
      double t3 = t2 * t;
      spline[i-1] = 0.5 * ((2.0*p1) + (-p0 + p2) * t + \
                    (2.0*p0 - 5.0*p1 + 4.0*p2 - p3) * t2 + \
                    (-p0 + 3.0*p1 - 3.0*p2 + p3) * t3);
    }
    return spline;
  }

}} // namespace scitbx::math

#endif
