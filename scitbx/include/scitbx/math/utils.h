#ifndef SCITBX_MATH_UTILS_H
#define SCITBX_MATH_UTILS_H

#include <cmath>

namespace scitbx { namespace math {

  inline
  int
  iround(double x)
  {
    if (x < 0) return static_cast<int>(x-0.5);
    return static_cast<int>(x+.5);
  }

  inline
  int
  iceil(double x) { return iround(std::ceil(x)); }

  inline
  int
  ifloor(double x) { return iround(std::floor(x)); }

}} // namespace scitbx::math

#endif // SCITBX_MATH_UTILS_H
