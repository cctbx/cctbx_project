#ifndef SCITBX_MATH_LINEAR_INTERPOLATION_H
#define SCITBX_MATH_LINEAR_INTERPOLATION_H

namespace scitbx { namespace math {

  template <typename FloatType>
  inline
  FloatType
  linear_interpolation(
    FloatType const& x,
    FloatType const& x1,
    FloatType const& x2,
    FloatType const& y1,
    FloatType const& y2)
  {
    FloatType dx = x2 - x1;
    FloatType dy = y2 - y1;
    return y1 + dy*(x-x1)/dx;
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_LINEAR_INTERPOLATION_H
