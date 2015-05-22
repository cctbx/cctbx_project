#ifndef SCITBX_FAST_APPROX_MATH_H
#define SCITBX_FAST_APPROX_MATH_H
#include <scitbx/constants.h>

namespace scitbx { namespace math {

template <typename FloatType>
FloatType
cos_table(
  af::const_ref<FloatType> const& table,
  FloatType arg,
  FloatType const& step,
  int const& n,
  bool interpolate)
{
  FloatType two_pi = scitbx::constants::two_pi;
  if(arg<0) arg = std::abs(arg);
  if(arg>two_pi)
    arg = arg - int(arg/two_pi)*two_pi;
  arg = arg/step;
  int k = scitbx::math::mod_positive(int(arg),n);
  if(interpolate) {
    FloatType y=table[k];
    return y+(table[scitbx::math::mod_positive(k+1,n)]-y)*(arg-k);
  }
  else {
    return table[k];
  }
}

template <typename FloatType>
FloatType sin_table(
         af::const_ref<FloatType> const& table,
         FloatType arg,
         FloatType const& step,
         int const& n,
         bool interpolate) {
  FloatType two_pi = scitbx::constants::two_pi;
  if(arg<0) {
    arg = std::abs(arg);
    if(arg>two_pi)
      arg = arg - int(arg/two_pi)*two_pi;
    arg = arg/step;
    int k = scitbx::math::mod_positive(int(arg),n);
    if(interpolate) {
      FloatType y=table[k];
      return -y-(table[scitbx::math::mod_positive(k+1,n)]-y)*(arg-k);
    }
    else {
      return -table[k];
    }
  }
  else {
    if(arg>two_pi)
      arg = arg - int(arg/two_pi)*two_pi;
    arg = arg/step;
    int k = scitbx::math::mod_positive(int(arg),n);
    if(interpolate) {
      FloatType y=table[k];
      return y+(table[scitbx::math::mod_positive(k+1,n)]-y)*(arg-k);
    }
    else {
      return table[k];
    }
  }
}

float approx_sqrt(float x) {
  unsigned int i = *(unsigned int*) &x;
  i  += 127 << 23;
  i >>= 1;
  return *(float*) &i;
}

}} // namespace scitbx::math

#endif // SCITBX_FAST_APPROX_MATH_H
