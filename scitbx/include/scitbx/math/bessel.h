#ifndef SCITBX_MATH_BESSEL_H
#define SCITBX_MATH_BESSEL_H

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

}}} // namespace scitbx::math::bessel

#endif // SCITBX_MATH_BESSEL_H
