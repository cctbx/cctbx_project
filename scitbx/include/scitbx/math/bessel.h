#ifndef SCITBX_MATH_BESSEL_H
#define SCITBX_MATH_BESSEL_H

#include <cmath>

namespace scitbx { namespace math {

//! Family of Bessel functions.
namespace bessel {

  //! Calculates the ratio I1(x)/I0(x).
  template <typename FloatType>
  FloatType
  i1_over_i0(FloatType const& x)
  {
    typedef FloatType f_t;
    f_t p[] = {1.0,3.5156229,3.0899424,1.2067292,0.2659732,0.360768E-1,
               0.45813E-2};
    f_t q[] = {0.39894228,0.1328592E-1,0.225319E-2,-0.157565E-2,0.916281E-2,
      -0.2057706E-1,0.2635537E-1,-0.1647633E-1,0.392377E-2};
    f_t pp[] = {0.5,0.87890594,0.51498869,0.15084934,0.2658733E-1,0.301532E-2,
        0.32411E-3};
    f_t qq[] = {0.39894228,-0.3988024E-1,-0.362018E-2,0.163801E-2,
      -0.1031555E-1, 0.2282967E-1,-0.2895312E-1,0.1787654E-1,-0.420059E-2};
    f_t be1 = 0;
    f_t be0 = 0;
    f_t abs_x = x;
    if (abs_x < 0) abs_x = -abs_x;
    if (abs_x < 3.75) {
      f_t y=x/3.75;
      y *= y;
      f_t pow_y_i = 1;
      for(int i=0; i<7; i++) {
        be0 += p[i]*pow_y_i;
        be1 += x*pp[i]*pow_y_i;
        pow_y_i *= y;
      }
    }
    else {
      f_t y=3.75/abs_x;
      f_t pow_y_i = 1;
      for(int i=0; i<9; i++) {
        be0 += q[i]*pow_y_i;
        be1 += qq[i]*pow_y_i;
        pow_y_i *= y;
      }
    }
    f_t result = be1/be0;
    if (x < 0.0 && result > 0.0) return -result;
    return result;
  }

  //! Calculates I0(x).
  template <typename FloatType>
  FloatType
  i0(FloatType const& x)
  {
    typedef FloatType f_t;
    f_t abs_x = x;
    if (abs_x < 0) abs_x = -abs_x;
    f_t bessel_i0;
    if (abs_x/3.75 < 1.0) {
       f_t y = x/3.75;
       y *= y;
       bessel_i0=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+
                 y*(0.2659732+y*(0.0360768+y*0.0045813)))));
    }
    else {
       f_t y = 3.75/abs_x;
       y=0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+
         y*(0.00916281+y*(-0.02057706+y*(0.02635537+
         y*(-0.01647633+y*0.00392377)))))));
       bessel_i0=y*std::exp(abs_x)/std::sqrt(abs_x);
    }
    return bessel_i0;
  }

  //! Calculates I1(x).
  template <typename FloatType>
  FloatType
  i1(FloatType const& x)
  {
    typedef FloatType f_t;
    f_t abs_x = x;
    if (abs_x < 0) abs_x = -abs_x;
    f_t ans;
    if (abs_x/3.75 < 1.0) {
      f_t y=x/3.75;
      y*=y;
      ans=abs_x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
          +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
    }
    else {
      f_t y=3.75/abs_x;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
          -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
          +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (std::exp(abs_x)/std::sqrt(abs_x));
    }
    f_t result = ans;
    if (x < 0.0 && result > 0.0) return -result;
    return result;
  }

  //! Calculates ln( I0(x) ).
  template <typename FloatType>
  FloatType
  ln_of_i0(FloatType const& x)
  {
    typedef FloatType f_t;
    f_t abs_x = x;
    if (abs_x < 0) abs_x = -abs_x;
    f_t bessel_lni0;
    if (abs_x/3.75 < 1.0) {
       f_t y = x/3.75;
       y *= y;
       bessel_lni0=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+
                 y*(0.2659732+y*(0.0360768+y*0.0045813)))));
       bessel_lni0 = std::log(bessel_lni0);
    }
    else {
       f_t y = 3.75/abs_x;
       y=0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+
         y*(0.00916281+y*(-0.02057706+y*(0.02635537+
         y*(-0.01647633+y*0.00392377)))))));
       bessel_lni0 = std::log(y) + abs_x - 0.5*std::log(abs_x);
    }
    return bessel_lni0;
  }

}}} // namespace scitbx::math::bessel

#endif // SCITBX_MATH_BESSEL_H
