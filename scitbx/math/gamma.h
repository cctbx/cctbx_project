// Peter Zwart, April 1st, 2005
#ifndef SCITBX_MATH_GAMMA_H
#define SCITBX_MATH_GAMMA_H

#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>

namespace scitbx { namespace math {

//! The Gamma function and some of his/her friends.
namespace gamma {

  //! Factorial/Gamma function.
  /*! Expression and coefficients obtained from:
      http://www.rskey.org/gamma.htm
   */
  template <typename FloatType>
  FloatType
  complete_lanczos(FloatType const& x)
  {
    if (x>=141.691) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "gamma::complete_lanczos(%.6g): domain error", x);
      throw error(buf);
    }
    typedef FloatType f_t;
    f_t q[] = {75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705,
               1168.92649479, 83.8676043424, 2.50662827511};
    f_t sum_part = q[0]*1.0; // sum_i=0^n=6 q[i]*x^n
    f_t prod_part = x;       // prod_i=0^n=6 (x+n)
    f_t z_tmp = 1.0;         // x^n
    for (int i=1; i<=6; i++){
      z_tmp *= x;
      sum_part += q[i]*z_tmp;
      prod_part *= (x+i);
    }
    return (sum_part/prod_part)
      *std::pow((x+5.5),(x+0.5))
      *std::exp(-x-5.5);
  }

  //! Minimax approach for arguments larger than 12.
  /*! http://www.netlib.no/netlib/specfun/gamma
   */
  template <typename FloatType>
  FloatType
  log_complete_minimax(FloatType const& x)
  {
    SCITBX_ASSERT(x > 12);
    typedef FloatType f_t;
    f_t sqrtpi = 0.9189385332046727417803297;
    f_t q[] = {-1.910444077728E-03,         8.4171387781295E-04,
               -5.952379913043012E-04,      7.93650793500350248E-04,
               -2.777777777777681622553E-03,8.333333333333333331554247E-02,
                5.7083835261E-03};
    f_t xsq = x*x;
    f_t sum = q[6];
    for (int i = 0;i<6;i++){
      sum = sum/xsq + q[i];
    }
    sum = sum/x -x + sqrtpi;
    sum = sum + (x-0.5)*std::log(x);
    return (sum);
  }

  template <typename FloatType>
  FloatType
  complete_minimax(FloatType const& x)
  {
    if (x>=171.624) {
      char buf[128];
      std::snprintf(buf, sizeof(buf), "gamma::complete_minimax(%.6g): domain error", x);
      throw error(buf);
    }
    return std::exp( log_complete_minimax(x) );
  }


  //! Gamma function with automatic choice of the best approximation.
  template <typename FloatType>
  FloatType
  complete(FloatType const& x, bool minimax=true)
  {
    if (minimax && x > 12) return complete_minimax(x);
    return complete_lanczos(x);
  }


  template <typename FloatType>
  FloatType
  log_complete(FloatType const& x, bool minimax=true)
  {
    if (minimax && x > 12) return log_complete_minimax(x);
    return std::log(complete_lanczos(x));
  }

  //! Incomplete gamma function using series expansion.
  template <typename FloatType>
  FloatType
  incomplete_series(FloatType const& a,
                    FloatType const& x,
                    unsigned max_iterations=500)
  {
    SCITBX_ASSERT(a > 0);
    SCITBX_ASSERT(x >= 0);
    if (x == 0) return 0;
    typedef FloatType f_t;
    f_t del = 1.0/a;
    f_t sum = 1.0/a;
    f_t eps = scitbx::math::floating_point_epsilon<FloatType>::get();
    for (unsigned i=1;i<=max_iterations;i++) {
      del *= x/(a+i);
      sum += del;
      if (std::fabs(del) < std::fabs(sum)*eps) {
        return sum*std::exp(-x+a*std::log(x)-std::log(complete(a)));
      }
    }
    char buf[256];
    std::snprintf(buf, sizeof(buf),
      "gamma::incomplete_series(a=%.6g, x=%.6g, max_iterations=%u)"
      " failed to converge", a, x, max_iterations);
    throw error(buf);
  }

  //! Incomplete gamma function using continued fraction method.
  template <typename FloatType>
  FloatType
  incomplete_continued_fraction(FloatType const& a,
                                FloatType const& x,
                                unsigned max_iterations=500)
  {
    SCITBX_ASSERT(a>0);
    SCITBX_ASSERT(x>=0);
    typedef FloatType f_t;
    static const f_t fpmin = 1.e-30;
    f_t eps = scitbx::math::floating_point_epsilon<FloatType>::get();
    f_t b = x+1.0-a;
    f_t c = 1.0/fpmin;
    f_t d = 1.0/b;
    f_t h = d;
    for (unsigned i=1;i<=max_iterations;i++){
      f_t an = -(i*(i-a));
      b += 2.0;
      d = an*d+b;
      if (std::fabs(d) < fpmin){ d = fpmin; }
      c = b+an/c;
      if (std::fabs(c) < fpmin){ c = fpmin; }
      d=1.0/d;
      f_t del = d*c;
      h *= del;
      if (std::fabs(del-1.0) < eps) {
        return 1.0-std::exp(-x+a*std::log(x)-std::log(complete(a)))*h;
      }
    }
    char buf[256];
    std::snprintf(buf, sizeof(buf),
      "gamma::incomplete_continued_fraction(a=%.6g, x=%.6g, max_iterations=%u)"
      " failed to converge", a, x, max_iterations);
    throw error(buf);
  }

  //! Normalized form of incomplete gamma function.
  /*! (1/\Gamma(a))int_{0}^{x}exp(-t) t^{a-1} dt
   */
  template <typename FloatType>
  FloatType
  incomplete(FloatType const& a,
             FloatType const& x,
             unsigned max_iterations=500)
  {
    SCITBX_ASSERT(a > 0);
    SCITBX_ASSERT(x >= 0);
    if (x < a+1.0) {
      return incomplete_series(a, x, max_iterations);
    }
    return incomplete_continued_fraction(a, x, max_iterations);
  }

  //! Complement of normalized form of incomplete gamma function.
  /*! 1-(1/\Gamma(a))int_{0}^{x}exp(-t) t^{a-1} dt
   */
  template <typename FloatType>
  FloatType
  incomplete_complement(FloatType const& a,
                        FloatType const& x,
                        unsigned max_iterations=500)
  {
    return 1 - incomplete(a, x, max_iterations);
  }

  /*! The exponential integral \integral_{z}^{\infty} \exp[-t] t^{-1} d t
   *  a modified form of AMS 55 5.1.53
   */
  template<typename FloatType>
  FloatType
  exponential_integral_e1z_lower_track(FloatType const& z)
  {
    SCITBX_ASSERT(z>0);
    SCITBX_ASSERT(z<=1);
    FloatType result=0;
    result = -std::log(z);
    FloatType a0,a1,a2,a3,a4,a5;
    a0=-0.57721566;
    a1= 0.99999193;
    a2=-0.24991055;
    a3= 0.05519968;
    a4=-0.00976004;
    a5= 0.00107857;
    result += a0 +
              a1*z +
              a2*z*z +
              a3*z*z*z +
              a4*z*z*z*z +
              a5*z*z*z*z*z;

    return(result);
  }

  /*! The exponential integral \integral_{z}^{\infty} \exp[-t] t^{-1} d t
   *  a modified form of AMS 55 5.1.56
   */
  template<typename FloatType>
  FloatType
  exponential_integral_e1z_upper_track(FloatType const& z)
  {
    SCITBX_ASSERT(z>=1);
    FloatType result=0;
    FloatType a1,a2,a3,a4;
    FloatType b1,b2,b3,b4;

    a1= 8.5733287401;
    a2=18.0590169730;
    a3= 8.6347608925;
    a4= 0.2677737343;

    b1= 9.5733223454;
    b2=25.6329561486;
    b3=21.0996530827;
    b4= 3.9584969228;

    FloatType top,bottom;
    top =   z*z*z*z +
         a1*z*z*z +
         a2*z*z +
         a3*z +
         a4;
    bottom =    z*z*z*z +
             b1*z*z*z +
             b2*z*z +
             b3*z +
             b4;
    result = std::log(top)-std::log(bottom);
    result = result-std::log(z) - z;
    result = std::exp( result );
    return( result );
  }

  /*! The exponential integral \integral_{z}^{\infty} \exp[-t] t^{-1} d t
   *  a modified form of AMS 55 5.1.56
   */
  template<typename FloatType>
  FloatType
  exponential_integral_e1z(FloatType const& z)
  {
     SCITBX_ASSERT(z>=0);
     if (z<1.0){
        return( exponential_integral_e1z_lower_track(z) );
     }
     return( exponential_integral_e1z_upper_track(z) );
  }

}}} // namespace scitbx::math::gamma

#endif // SCITBX_MATH_GAMMA_H
