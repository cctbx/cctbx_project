#ifndef SCITBX_MATH_BESSEL_H
#define SCITBX_MATH_BESSEL_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <cmath>

#if !(defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 2)
# include <boost/math/special_functions/bessel.hpp>
# define SCITBX_MATH_BESSEL_HAS_SPHERICAL
#endif

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

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  i1_over_i0( scitbx::af::const_ref<FloatType> const& x )
  {
    SCITBX_ASSERT( x.size()>0 );
    scitbx::af::shared<FloatType> result;
    for (int ii=0;ii<x.size();ii++){
      result.push_back( i1_over_i0(x[ii]) );
    }
    return(result);
  }


  //! Calculates the inverse function of the ratio I1(x)/I0(x).
  /*! Implementation based on clipper::Util::invsim() by Kevin Cowtan.
      Warning: for x > 10 the results are not very accurate.
   */
  template <typename FloatType>
  FloatType
  inverse_i1_over_i0(const FloatType& x)
  {
    typedef FloatType f_t;
    f_t x0 = std::fabs(x);
    f_t a0 = -7.107935*x0;
    f_t a1 = 3.553967-3.524142*x0;
    f_t a2 = 1.639294-2.228716*x0;
    f_t a3 = 1.0-x0;
    f_t w = a2/(3.0*a3);
    f_t p = a1/(3.0*a3)-w*w;
    f_t q = -w*w*w+0.5*(a1*w-a0)/a3;
    f_t d = std::sqrt(q*q+p*p*p);
    f_t q1 = q + d;
    f_t q2 = q - d;
    f_t r1 = std::pow(std::fabs(q1), 1.0/3.0);
    f_t r2 = std::pow(std::fabs(q2), 1.0/3.0);
    if (x >= 0.0) return  (((q1>0.0)? r1 : -r1) + ((q2>0.0)? r2 : -r2) - w);
                  return -(((q1>0.0)? r1 : -r1) + ((q2>0.0)? r2 : -r2) - w);
  }

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  inverse_i1_over_i0( scitbx::af::const_ref<FloatType> const& x )
  {
    SCITBX_ASSERT( x.size()>0 );
    scitbx::af::shared<FloatType> result;
    for (int ii=0;ii<x.size();ii++){
      result.push_back( inverse_i1_over_i0(x[ii]) );
    }
    return(result);
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

  template <typename FloatType>
  FloatType
  ei1(FloatType const& x)
  {
    // a quick and dirty approximation
    FloatType t;
    t = x/(1.0+x);
    /*


    Final set of parameters            Asymptotic Standard Error
    =======================            ==========================

    a               = 0.47735          +/- 0.0003554    (0.07446%)
    b               = 0.175945         +/- 0.004666     (2.652%)
    c               = -2.82479         +/- 0.02347      (0.831%)
    d               = 8.74696          +/- 0.05844      (0.6681%)
    e               = -14.2023         +/- 0.07694      (0.5417%)
    f               = 10.9149          +/- 0.05134      (0.4704%)
    g               = -3.14135         +/- 0.01368      (0.4356%)


    correlation matrix of the fit parameters:

                    a      b      c      d      e      f      g
    a               1.000
    b              -0.974  1.000
    c               0.930 -0.988  1.000
    d              -0.887  0.965 -0.993  1.000
    e               0.847 -0.939  0.979 -0.996  1.000
    f              -0.812  0.912 -0.962  0.986 -0.997  1.000
    g               0.781 -0.887  0.943 -0.974  0.990 -0.998  1.000
    gnuplot> f(x) = sqrt(1-x)*exp(x)*(x*(a+x*(b+x*(c+x*(d+x*(e+x*(f+g*x)))))))
    */
    FloatType a=  0.47735;
    FloatType b=  0.175945;
    FloatType c= -2.82479;
    FloatType d=  8.74696;
    FloatType e=-14.2023;
    FloatType f= 10.9149;
    FloatType g= -3.14135;

    FloatType result;
    result = std::sqrt(1-t)*std::exp(t)*(t*(a+t*(b+t*(c+t*(d+t*(e+t*(f+g*t)))))));
    return result;
  }


  template <typename FloatType>
  FloatType
  ei0(FloatType const& x)
  {
    // A quick and dirty approximation
    FloatType t;
    t = x/(1.0+x);
    /*
      degrees of freedom (ndf) : 39
      rms of residuals      (stdfit) = sqrt(WSSR/ndf)      : 0.000162321
      variance of residuals (reduced chisquare) = WSSR/ndf : 2.6348e-08

      Final set of parameters            Asymptotic Standard Error
      =======================            ==========================

      b               = -1.51857         +/- 0.002712     (0.1786%)
      c               = 0.862203         +/- 0.01986      (2.304%)
      d               = -1.11554         +/- 0.05248      (4.705%)
      e               = 1.72229          +/- 0.05872      (3.41%)
      f               = -0.804154        +/- 0.02353      (2.926%)


      correlation matrix of the fit parameters:

                      b      c      d      e      f
      b               1.000
      c              -0.973  1.000
      d               0.922 -0.986  1.000
      e              -0.868  0.956 -0.992  1.000
      f               0.819 -0.923  0.973 -0.995  1.000
      gnuplot> plot 'pp' using 2:3, f(x)
      gnuplot>  f(x) = sqrt(1-x)*exp(x)*(1+x*(b+x*(c+x*(d+x*(e+f*x)))))
    */

    FloatType a=  1.0;
    FloatType b= -1.51857;
    FloatType c=  0.862203;
    FloatType d= -1.11554;
    FloatType e=  1.72229;
    FloatType f= -0.804154;

    FloatType result;
    result = std::sqrt(1-t)*std::exp(t)*(a+t*(b+t*(c+t*(d+t*(e+f*t)))))  ;
    return result;
  }

#if defined(SCITBX_MATH_BESSEL_HAS_SPHERICAL)
  template <typename FloatType>
  FloatType
  spherical_bessel(int const& l, FloatType const& x)
  {
    using boost::math::sph_bessel;
    return( sph_bessel(l,x) );
  }

  template <typename FloatType>
  scitbx::af::shared< FloatType >
  spherical_bessel_array(int const& l, scitbx::af::shared< FloatType >const& x)
  {
    using boost::math::sph_bessel;
    scitbx::af::shared< FloatType > result;
    for (int ii=0;ii<x.size();ii++){
      result.push_back( spherical_bessel( l, x[ii] ) );
    }
    return( result );
  }
#endif

}}} // namespace scitbx::math::bessel

#endif // SCITBX_MATH_BESSEL_H
