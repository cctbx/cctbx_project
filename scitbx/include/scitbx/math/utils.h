#ifndef SCITBX_MATH_UTILS_H
#define SCITBX_MATH_UTILS_H

#include <cmath>

namespace scitbx { namespace math {

  template <typename FloatType,
            typename SignedIntType>
  struct float_int_conversions
  {
    static
    inline
    SignedIntType
    iround(FloatType const& x)
    {
      if (x < 0) return static_cast<SignedIntType>(x-0.5);
      return static_cast<SignedIntType>(x+.5);
    }

    static
    inline
    SignedIntType
    iceil(FloatType const& x) { return iround(std::ceil(x)); }

    static
    inline
    SignedIntType
    ifloor(FloatType const& x) { return iround(std::floor(x)); }
  };

  inline
  int
  iround(double x) { return float_int_conversions<double,int>::iround(x); }

  inline
  int
  iceil(double x) { return float_int_conversions<double,int>::iceil(x); }

  inline
  int
  ifloor(double x) { return float_int_conversions<double,int>::ifloor(x); }

}} // namespace scitbx::math

#endif // SCITBX_MATH_UTILS_H
