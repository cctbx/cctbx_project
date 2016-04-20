#ifndef SCITBX_MATH_LINEAR_INTERPOLATION_H
#define SCITBX_MATH_LINEAR_INTERPOLATION_H

#include <scitbx/error.h>

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
  };

  template <typename FloatType>
  inline
  FloatType
  linear_interpolation_2d(
    FloatType const& x1,
    FloatType const& y1,
    FloatType const& x2,
    FloatType const& y2,
    FloatType const& v1,
    FloatType const& v2,
    FloatType const& v3,
    FloatType const& v4,
    FloatType const& xx,
    FloatType const& yy)
  {
    // Perform interpolation for 2d table.
    // Y ^   v3       v32
    //   |    *-------|------------* x2,y2,v2
    //   |    |                    |
    //   |    |                    |
    //   |    |                    |
    //   |    |       *            |
    //   |    |     xx,yy          |
    //   |    |                    |
    //   |    |                    |
    //   |    |                    |
    //   |    *-------|------------*
    //   | x1,y1,v1   v14           v4
    //   |
    //  -|---------------------------> X
    //
    // https://en.wikipedia.org/wiki/Bilinear_interpolation
    //  Returns interpolated value at xx,yy
    //
    SCITBX_ASSERT(x1 < x2);
    SCITBX_ASSERT(y1 < y2);
    SCITBX_ASSERT(x1 <= xx);
    SCITBX_ASSERT(xx <= x2);
    SCITBX_ASSERT(y1 <= yy);
    SCITBX_ASSERT(yy <= y2);
    // first we do interpolation in X direction
    FloatType v14 = linear_interpolation(xx,x1,x2,v1,v4);
    FloatType v32 = linear_interpolation(xx,x1,x2,v3,v2);
    // then we do interpolation in Y direction
    FloatType result = linear_interpolation(yy,y1,y2,v14,v32);
    return result;
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_LINEAR_INTERPOLATION_H
