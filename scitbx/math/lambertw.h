#ifndef SCITBX_MATH_LAMBERTW_H
#define SCITBX_MATH_LAMBERTW_H

#include <scitbx/math/floating_point_epsilon.h>

namespace scitbx { namespace math {

  //! Lambert's W function: w(x)*e^(w(x))=x
  /*! Based on the implementation by Gunter Kuhnle, gk@uni-leipzig.de

      Algorithm originally developed by
      Keith Briggs, Department of Plant Sciences,
      e-mail:kmb28@cam.ac.uk

      http://more.btexact.com/people/briggsk2/W-ology.html

      See also:
        http://documents.wolfram.com/mathematica/functions/ProductLog
   */
  template <typename FloatType>
  FloatType
  lambertw(
    FloatType const& x,
    unsigned max_iterations=100)
  {
    static FloatType eps = 0;
    if (eps == 0) {
      eps = scitbx::math::floating_point_epsilon<FloatType>::get();
    }
    if (x < -std::exp(static_cast<FloatType>(-1))) {
      throw std::runtime_error("lambertw(x) domain error: x < -exp(-1)");
    }
    if (std::fabs(x) <= eps) {
      return x;
    }
    FloatType w;
    if (x < 1.0) {
      FloatType p = std::sqrt(2.0 * (std::exp(1.0) * x + 1.0));
      FloatType pp = p * p;
      w = -1.0 + p - pp / 3.0 + 11.0 / 72.0 * pp * p;
    }
    else {
      w = std::log(x);
    }
    if (x > 3) {
      SCITBX_ASSERT(w>0);
      w -= std::log(w);
    }
    for(unsigned i=0;i<max_iterations;i++) {
      FloatType e = std::exp(w);
      FloatType t = w * e - x;
      t /= (e * (w + 1.0) - 0.5 * (w + 2.0) * t / (w + 1.0));
      w -= t;
      if (std::fabs(t) < eps * (1.0 + std::fabs(w))) {
        return w;
      }
    }
    throw std::runtime_error("lambertw error: iteration did not converge");
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_LAMBERTW_H
