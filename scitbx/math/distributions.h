#ifndef SCITBX_MATH_DISTRIBUTIONS_H
#define SCITBX_MATH_DISTRIBUTIONS_H

#include <scitbx/array_family/shared.h>
#include <boost/math/distributions.hpp>

/*! Extension of boost::math statistical distributions.
    See also:
    http://www.boost.org/libs/math/doc/sf_and_dist/html/math_toolkit/dist.html
 */
namespace scitbx { namespace math {

  /*! An array of quantiles of size n.
      For a given size n, then for k = 1,...,n,
      the quantiles chosen are (k-.5)/n.
   */
  template <typename FloatType, class Distribution>
  af::shared<FloatType>
  quantiles(const Distribution& dist, std::size_t n)
  {
    af::shared<FloatType> result(n);
    for (std::size_t i=0;i<n;i++) {
      result[i] = boost::math::quantile(
        dist, (i+.5)/n);
    }
    return result;
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_DISTRIBUTIONS_H
