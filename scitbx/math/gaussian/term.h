#ifndef SCITBX_MATH_GAUSSIAN_TERM_H
#define SCITBX_MATH_GAUSSIAN_TERM_H

#include <scitbx/math/erf.h>
#include <scitbx/constants.h>
#include <cmath>

namespace scitbx { namespace math { namespace gaussian {

  //! Gaussian term: a Exp[-b x^2]
  template <typename FloatType=double>
  struct term
  {
    //! Default constructor. Some data members are not initialized!
    term() {}

    //! Definition of coefficients a and b.
    term(FloatType const& a_, FloatType const& b_)
    :
      a(a_),
      b(b_)
    {}

    //! Gaussian at the point x, given x^2.
    FloatType
    at_x_sq(FloatType const& x_sq) const
    {
      return a * std::exp(-b * x_sq);
    }

    //! Gaussian at the point x.
    FloatType
    at_x(FloatType const& x) const
    {
      return at_x_sq(x * x);
    }

    //! Analytical gradient w.r.t. x at the point x.
    FloatType
    gradient_dx_at_x(FloatType const& x) const
    {
      return -2*a*b*x/std::exp(b*x*x);
    }

    //! Analytical integral dx from 0 to the point x.
    FloatType
    integral_dx_at_x(
      FloatType const& x,
      FloatType const& b_min_for_erf_based_algorithm=1e-3)
    {
      using scitbx::math::erf;
      static const FloatType sqrt_pi = std::sqrt(scitbx::constants::pi);
      if (b == 0) return a * x;
      if (b > b_min_for_erf_based_algorithm) {
        /* Mathematica:
             f = a Exp[-b x^2]
             Integrate[f,x]
         */
        FloatType sqrt_b = std::sqrt(b);
        return a*sqrt_pi*erf(sqrt_b*x)/(2*sqrt_b);
      }
      /* Mathematica:
           f = a Exp[-b x^2]
           Series[Integrate[f,x], {x,0,20}]
         Formula for the denominator of the series expansion: (2n+1)*n!
         Encyclopedia of Integer Sequences ID Number: A007680
       */
      FloatType bxx = b * x * x;
      FloatType part = 1;
      FloatType result = 1;
      FloatType prev_result = result;
      unsigned n = 0;
      unsigned tnp1 = 1;
      while (true) {
        n++;
        tnp1 += 2;
        part *= bxx / n;
        result -= part / tnp1;
        if (result == prev_result) break;
        prev_result = result;
        n++;
        tnp1 += 2;
        part *= bxx / n;
        result += part / tnp1;
        if (result == prev_result) break;
        prev_result = result;
      }
      return a * x * result;
    }

    //! Analytical gradients w.r.t. a and b at the point x, given x^2.
    term
    gradients_d_ab_at_x_sq(FloatType const& x_sq) const
    {
      FloatType gr_a = std::exp(-b * x_sq);
      return term(gr_a, -a * x_sq * gr_a);
    }

    FloatType a;
    FloatType b;
  };

}}} // scitbx::math::gaussian

#endif // SCITBX_MATH_GAUSSIAN_TERM_H
