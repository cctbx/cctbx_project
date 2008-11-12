#ifndef SCITBX_MATH_MOD_H
#define SCITBX_MATH_MOD_H

#include <cmath>

namespace scitbx { namespace math {

  template <typename IntType>
  inline
  IntType
  mod_positive(IntType ix, IntType const& iy)
  {
    if (iy > 0) {
      ix %= iy;
      if (ix < 0) ix += iy;
    }
    return ix;
  }

  template <typename IntType>
  inline
  IntType
  mod_short(IntType ix, IntType const& iy)
  {
        ix = mod_positive(ix, iy);
    if (ix > iy / 2)
        ix -= iy;
    return ix;
  }

  template <typename FloatType>
  inline FloatType
  fmod_short(FloatType const& x, FloatType const& y)
  {
    FloatType result = std::fmod(x, y);
    if (result < 0) result += y;
    if (result > y/2) result -= y;
    return result;
  }

}} // namespace scitbx::math

#endif // GUARD
