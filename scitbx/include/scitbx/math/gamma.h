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
      std::sprintf(buf, "gamma::complete_lanczos(%.6g): domain error", x);
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
  complete_minimax(FloatType const& x)
  {
    if (x>=171.624) {
      char buf[128];
      std::sprintf(buf, "gamma::complete_minimax(%.6g): domain error", x);
      throw error(buf);
    }
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
    return std::exp(sum);
  }

  //! Gamma function with automatic choice of the best approximation.
  template <typename FloatType>
  FloatType
  complete(FloatType const& x, bool minimax=true)
  {
    if (minimax && x > 12) return complete_minimax(x);
    return complete_lanczos(x);
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
    std::sprintf(buf,
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
    std::sprintf(buf,
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

}}} // namespace scitbx::math::gamma

#endif // SCITBX_MATH_GAMMA_H
