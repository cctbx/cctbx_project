// Peter Zwart, April 1st, 2005
#ifndef SCITBX_MATH_GAMMA_H
#define SCITBX_MATH_GAMMA_H

#include <iostream>
#include <cmath>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/constants.h>


namespace scitbx { namespace math {

//! The Gamma function and some of his/her friends

  static const double fpmin = 1E-30;
  static const double xbig = 171.624;
  static const double xinf = 1.97E307;

  /*
   * Calculates the factorial/Gamma function
   * Expression and coefficients obtained from
   * http://www.rskey.org/gamma.htm
   */

  template <typename FloatType>
  FloatType
  complete_gamma_lanczos(FloatType const& x)
  {

    typedef FloatType f_t;

    f_t q[] = {75122.6331530, 80916.6278952, 36308.2951477, 8687.24529705,
               1168.92649479, 83.8676043424, 2.50662827511};

    f_t sum_part = 0.0;  // sum_i=0^n=6 q[i]*x^n
    f_t prod_part = 0.0; // prod_i=0^n=6 (x+n)
    f_t z_tmp = 1.0;     // x^n

    sum_part = q[0]*1.0;
    prod_part = x+0.0;

    for (int i=1; i<=6; i++){
      z_tmp = z_tmp*x;
      sum_part = sum_part + q[i]*z_tmp;
      prod_part = prod_part*(x+i);
    }
    f_t result;
    result = (sum_part/prod_part)
      *std::pow((x+5.5),(x+0.5))
      *std::exp(-x-5.5);
    return result;
  }

  /*
   * minimax approach for arguments larger then 12.
   * http://www.netlib.no/netlib/specfun/gamma
   */
  template <typename FloatType>
  FloatType
  complete_gamma_minimax(FloatType const& x)
  {

    SCITBX_ASSERT (x>12);

    typedef FloatType f_t;
    f_t sqrtpi = 0.9189385332046727417803297;
    f_t q[] = {-1.910444077728E-03,         8.4171387781295E-04,
               -5.952379913043012E-04,      7.93650793500350248E-04,
               -2.777777777777681622553E-03,8.333333333333333331554247E-02,
                5.7083835261E-03};
    f_t result = 0;
    if (x<xbig){
      f_t xsq = x*x;
      f_t sum = q[6];
      for (int i = 0;i<6;i++){
        sum = sum/xsq + q[i];
      }
      sum = sum/x -x + sqrtpi;
      sum = sum + (x-0.5)*std::log(x);
      result = std::exp(sum);
      return(result);
    }
    if (x>=xbig){
      return(xinf);
    }
    // never end up here
    return(EXIT_FAILURE);
  }

  /*
   * Compute the gamma function with the ability to
   * choose which approximation is needed.
   */

  template <typename FloatType>
  FloatType
  complete_gamma(FloatType const& x, bool minimax=true)
  {
    typedef FloatType f_t;
    if (x>12) {
      if (minimax){
        return ( complete_gamma_minimax(x) );
      }
      if (!minimax){
        return ( complete_gamma_lanczos(x) );
      }
    } else {
      return ( complete_gamma_lanczos(x) );
    }

    // never end up here
    return(EXIT_FAILURE);
  }

    /*
     * Computes incomplete gamma function using series expansion.
     */

  template <typename FloatType>
  FloatType
  incomplete_gamma_series(FloatType const& a,
                          FloatType const& x,
                          unsigned max_count=500)
  {

    typedef FloatType f_t;
    f_t eps = scitbx::math::floating_point_epsilon<FloatType>::get();
    f_t result=0.0;


    SCITBX_ASSERT(a>0);
    SCITBX_ASSERT(x>=0);

    if (x>0){
      f_t ap = a;
      f_t del = 1.0/a;
      f_t sum = 1.0/a;

      for (int i=1;i<max_count;i++){
        ++ap;
        del *= x/ap;
        sum += del;
        if (std::fabs(del) < std::fabs(sum)*eps){
          result = sum*std::exp(-x+a*std::log(x)-std::log(complete_gamma(a)));
          return(result);
        }

      }
      std::cout << "Incomplete gamma series failed to converge " << std::endl;
      std::cout << "Returning best guess " << std::endl;
      return(result);

    }

    return(0);
  }


    /*
     * Computes the incomplete gamma function using continued fraction method
     */


  template <typename FloatType>
  FloatType
  incomplete_gamma_continued_fraction(FloatType const& a,
                                      FloatType const& x,
                                      unsigned max_iterations=500)
  {

    SCITBX_ASSERT(a>0);
    SCITBX_ASSERT(x>=0);

    typedef FloatType f_t;
    f_t eps = scitbx::math::floating_point_epsilon<FloatType>::get();
    f_t an,b,c,d,del,h;
    f_t result;
    int i;
    b = x+1.0-a;
    c = 1.0/fpmin;
    d = 1.0/b;
    h = d;
    for (i=1;i<max_iterations;i++){
      an = -i*(i-a);
      b +=2.0;
      d = an*d+b;
      if (std::fabs(d) < fpmin){ d = fpmin; }
      c = b+an/c;
      if (std::fabs(c) < fpmin){ c = fpmin; }
      d=1.0/d;
      del = d*c;
      h *= del;
      if (std::fabs(del-1.0) <eps){ break; }
    }
    if (i>=max_iterations) {
      std::cout << "Continued fractions for incomplete gamma function failed to converge" << std::endl;
      std::cout << "Returning best guess" << std::endl;
    }

    result = 1.0-std::exp(-x+a*std::log(x)-std::log(complete_gamma(a)))*h;
    return(result);
  }


  /* Computes the incomplete gamma function in it's normalised form:
   * (1/\Gamma(a))int_{0}^{x}exp(-t) t^{a-1} dt
   */

  template <typename FloatType>
  FloatType
  incomplete_gamma(FloatType const& a,
                   FloatType const& x,
                   unsigned max_iterations=500)
  {
    SCITBX_ASSERT (x >= 0.0);
    SCITBX_ASSERT (a >  0.0);
    if (x < a+1.0) {// Use series expansion.
      return incomplete_gamma_series(a,x,max_iterations);
    } else { // Use cont. fraction.
      return incomplete_gamma_continued_fraction(a,x);
    }
    // never end up here
    return(EXIT_FAILURE);
  }

  /* Computes the complement of the incomplete gamma function in it's normalised form:
   * 1-(1/\Gamma(a))int_{0}^{x}exp(-t) t^{a-1} dt
   */

  template <typename FloatType>
  FloatType
  incomplete_gamma_complement(FloatType const& a,
                              FloatType const& x,
                              unsigned max_iterations=500)
  {
    SCITBX_ASSERT(x >= 0.0);
    SCITBX_ASSERT(a >  0.0);
    return (1.0 - incomplete_gamma(a,x,max_iterations));
  }



}} // namespace scitbx::math::gamma

#endif // SCITBX_MATH_GAMMA_H
