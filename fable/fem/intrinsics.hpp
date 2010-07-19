#ifndef FEM_INTRINSICS_HPP
#define FEM_INTRINSICS_HPP

#include <fem/str_ref.hpp>
#include <cmath>
#include <math.h> // C99 lround not in std namespace

namespace fem {

  template <typename T>
  inline
  int
  if_arithmetic(
    T const& value)
  {
    if (value == 0) return 0;
    if (value > 0) return 1;
    return -1;
  }

  template <typename V, typename S>
  inline
  V
  sign(
    V const& value,
    S const& sign_source)
  {
    if (sign_source < 0) {
      if (value > 0) return -value;
    }
    else if (value < 0) {
      return -value;
    }
    return value;
  }

  template <typename T>
  inline
  float
  real(
    T const& value) { return static_cast<float>(value); }

  template <typename T>
  inline
  double
  dble(
    T const& value) { return static_cast<double>(value); }

  inline
  float
  sqrt(
    float const& x)
  {
    return std::sqrt(x);
  }

  inline
  double
  sqrt(
    double const& x)
  {
    return std::sqrt(x);
  }

  inline
  double
  dsqrt(
    double const& x)
  {
    return std::sqrt(x);
  }

  inline
  float
  cos(
    float const& angle)
  {
    return std::cos(angle);
  }

  inline
  double
  cos(
    double const& angle)
  {
    return std::cos(angle);
  }

  inline
  float
  sin(
    float const& angle)
  {
    return std::sin(angle);
  }

  inline
  double
  sin(
    double const& angle)
  {
    return std::sin(angle);
  }

  inline
  float
  exp(
    float const& x)
  {
    return std::exp(x);
  }

  inline
  double
  exp(
    double const& x)
  {
    return std::exp(x);
  }

  inline
  double
  dexp(
    double const& x)
  {
    return std::exp(x);
  }

  inline
  float
  alog10(
    float const& x)
  {
    return std::log10(x);
  }

  inline
  double
  alog10(
    double const& x)
  {
    return std::log10(x);
  }

  template <typename T>
  inline
  int
  fint(
    T const& val)
  {
    return static_cast<int>(val);
  }

  template <typename T>
  inline
  int
  aint(
    T const& val)
  {
    return static_cast<int>(val);
  }

  template <typename T>
  inline
  float
  ffloat(
    T const& val)
  {
    return static_cast<float>(val);
  }

  inline
  int
  mod(
    int const& v1,
    int const& v2) { return v1 % v2; }

  inline
  float
  amod(
    float const& v1,
    float const& v2) { return std::fmod(v1, v2); }

  inline
  float
  mod(
    float const& v1,
    float const& v2) { return std::fmod(v1, v2); }

  inline
  double
  dmod(
    double const& v1,
    double const& v2) { return std::fmod(v1, v2); }

  inline
  double
  mod(
    double const& v1,
    double const& v2) { return std::fmod(v1, v2); }

  inline
  int
  iabs(
    int const& v) { return std::abs(v); }

  inline
  double
  dabs(
    double const& v) { return std::abs(v); }

  inline
  int
  min0(
    int const& v1,
    int const& v2) { return std::min(v1, v2); }

  inline
  int
  min(
    int const& v1,
    int const& v2) { return std::min(v1, v2); }

  inline
  float
  amin1(
    float const& v1,
    float const& v2) { return std::min(v1, v2); }

  inline
  float
  min(
    float const& v1,
    float const& v2) { return std::min(v1, v2); }

  inline
  double
  dmin1(
    double const& v1,
    double const& v2) { return std::min(v1, v2); }

  inline
  double
  min(
    double const& v1,
    double const& v2) { return std::min(v1, v2); }

  inline
  float
  amin0(
    int const& v1,
    float const& v2) { return std::min(static_cast<float>(v1), v2); }

  inline
  float
  min(
    int const& v1,
    float const& v2) { return std::min(static_cast<float>(v1), v2); }

  inline
  float
  min1(
    float const& v1,
    int const& v2) { return std::min(v1, static_cast<float>(v2)); }

  inline
  float
  min(
    float const& v1,
    int const& v2) { return std::min(v1, static_cast<float>(v2)); }

  inline
  int
  min(
    int const& v1,
    int const& v2,
    int const& v3)
  {
    return min(v1, min(v2, v3));
  }

  inline
  int
  min(
    int const& v1,
    int const& v2,
    int const& v3,
    int const& v4)
  {
    return min(v1, min(v2, v3, v4));
  }

  inline
  float
  min(
    float const& v1,
    float const& v2,
    float const& v3)
  {
    return min(v1, min(v2, v3));
  }

  inline
  float
  min(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4)
  {
    return min(v1, min(v2, v3, v4));
  }

  inline
  double
  min(
    double const& v1,
    double const& v2,
    double const& v3)
  {
    return min(v1, min(v2, v3));
  }

  inline
  double
  min(
    double const& v1,
    double const& v2,
    double const& v3,
    double const& v4)
  {
    return min(v1, min(v2, v3, v4));
  }

  inline
  int
  max0(
    int const& v1,
    int const& v2) { return std::max(v1, v2); }

  inline
  int
  max(
    int const& v1,
    int const& v2) { return std::max(v1, v2); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2) { return std::max(v1, v2); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2,
    float const& v3) { return amax1(amax1(v1, v2), v3); }

  inline
  float
  amax1(
    float const& v1,
    float const& v2,
    float const& v3,
    float const& v4) { return amax1(amax1(v1, v2, v3), v4); }

  inline
  float
  max(
    float const& v1,
    float const& v2) { return std::max(v1, v2); }

  inline
  double
  dmax1(
    double const& v1,
    double const& v2) { return std::max(v1, v2); }

  inline
  double
  max(
    double const& v1,
    double const& v2) { return std::max(v1, v2); }

  inline
  float
  amax0(
    int const& v1,
    float const& v2) { return std::max(static_cast<float>(v1), v2); }

  inline
  float
  max(
    int const& v1,
    float const& v2) { return std::max(static_cast<float>(v1), v2); }

  inline
  float
  max1(
    float const& v1,
    int const& v2) { return std::max(v1, static_cast<float>(v2)); }

  inline
  float
  max(
    float const& v1,
    int const& v2) { return std::max(v1, static_cast<float>(v2)); }

  inline
  double
  max(
    double const& v1,
    double const& v2,
    double const& v3) { return max(max(v1, v2), v3); }

  template <typename T>
  inline
  T
  pow2(
    T const& base) { return base * base; }

  template <typename T>
  inline
  T
  pow3(
    T const& base) { return base * base * base; }

  template <typename T>
  inline
  T
  pow4(
    T const& base) { T p2 = pow2(base); return p2 * p2; }

  inline
  int
  pow(
    int const& base,
    int const& exponent)
  {
    if (exponent < 0) return 0;
    int result = 1;
    for(int i=0;i<exponent;i++) {
      result *= base;
    }
    return result;
  }

  inline
  float
  pow(
    int const& base,
    float const& exponent)
  {
    return std::pow(static_cast<float>(base), exponent);
  }

  using std::abs;
  using std::log;
  using std::pow;
  using std::acos;
  using std::atan2;

  inline
  int
  ichar(
    str_cref c)
  {
    if (c.len() == 0) {
      std::ostringstream o;
      o << "ichar() argument must be a one-character string,"
        << " but actual string length is " << c.len() << ".";
      throw std::runtime_error(o.str());
    }
    return static_cast<int>(c[0]);
  }

  inline
  str<1>
  fchar(
    int i) { return str<1>(static_cast<char>(i)); }

  inline
  int
  len_trim(
    str_cref c)
  {
    return static_cast<int>(
      utils::find_trailing_blank_padding(c.elems(), c.len()));
  }

  inline
  int
  nint(
    float const& value) { return static_cast<int>(lroundf(value)); }

  inline
  int
  nint(
    double const& value) { return static_cast<int>(lround(value)); }

} // namespace fem

#endif // GUARD
