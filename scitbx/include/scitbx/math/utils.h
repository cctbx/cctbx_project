#ifndef SCITBX_MATH_UTILS_H
#define SCITBX_MATH_UTILS_H

namespace scitbx { namespace math {

  inline
  int
  iround(double x)
  {
    if (x < 0) return static_cast<int>(x-0.5);
    return static_cast<int>(x+.5);
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_UTILS_H
