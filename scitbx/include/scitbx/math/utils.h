#ifndef SCITBX_MATH_UTILS_H
#define SCITBX_MATH_UTILS_H

#include <cmath>

namespace scitbx { namespace math {

  template <typename NumType1, typename NumType2>
  inline
  NumType1&
  update_min(NumType1& m, NumType2 const& x)
  {
    if (m > x) m = x;
    return m;
  }

  template <typename NumType1, typename NumType2>
  inline
  NumType1&
  update_max(NumType1& m, NumType2 const& x)
  {
    if (m < x) m = x;
    return m;
  }

  template <typename FloatType>
  FloatType
  round(FloatType x, int n_digits=0)
  {
    // based on Python/bltinmodule.c: builtin_round()
    FloatType f = 1;
    int i = n_digits;
    if (i < 0) i = -i;
    while (--i >= 0) f *= 10;
    if (n_digits < 0) x /= f;
    else              x *= f;
    if (x >= 0) x = std::floor(x + 0.5);
    else        x = std::ceil(x - 0.5);
    if (n_digits < 0) x *= f;
    else              x /= f;
    return x;
  }

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

    static
    inline
    SignedIntType
    nearest_integer(FloatType const& x)
    {
      SignedIntType i = static_cast<SignedIntType>(x);
      if (x >= 0) {
        if (x - static_cast<FloatType>(i) >= .5) i++;
      }
      else {
        if (x - static_cast<FloatType>(i) <= -.5) i--;
      }
      return i;
    }
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

  inline
  int
  nearest_integer(double x)
  {
    return float_int_conversions<double,int>::nearest_integer(x);
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_UTILS_H
