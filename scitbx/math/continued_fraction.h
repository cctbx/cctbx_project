#ifndef SCITBX_MATH_CONTINUED_FRACTION
#define SCITBX_MATH_CONTINUED_FRACTION

#include <boost/rational.hpp>
#include <scitbx/math/utils.h>

namespace scitbx { namespace math {

template <typename IntType>
class continued_fraction
{
public:
  typedef IntType integral_type;
  typedef boost::rational<IntType> rational_type;

  continued_fraction(integral_type a0)
    : p_i_minus_1(1), q_i_minus_1(0), p_i(a0), q_i(1)
  {}

  void append(integral_type a) {
    integral_type p_i_plus_1 = a*p_i + p_i_minus_1;
    integral_type q_i_plus_1 = a*q_i + q_i_minus_1;
    p_i_minus_1 = p_i; p_i = p_i_plus_1;
    q_i_minus_1 = q_i; q_i = q_i_plus_1;
  }

  rational_type operator()() {
    return rational_type(p_i, q_i);
  }

  template<typename FloatType>
  FloatType as_real() {
    FloatType num = p_i, den = q_i;
    return num/den;
  }

  template <typename FloatType>
  static continued_fraction
  from_real(
    FloatType x, FloatType eps=std::numeric_limits<FloatType>::epsilon())
  {
    integral_type a = ifloor(x);
    continued_fraction result(a);
    FloatType y = x;
    while (std::abs(x - result.as_real<FloatType>()) > eps) {
      FloatType r = y-a;
      y = 1/r;
      a = ifloor(y);
      result.append(a);
    }
    return result;
  }

private:
  integral_type p_i, q_i, p_i_minus_1, q_i_minus_1;
};

}}

#endif // GUARD
