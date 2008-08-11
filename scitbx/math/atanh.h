#ifndef SCITBX_MATH_ATANH_H
#define SCITBX_MATH_ATANH_H

#if defined(__GNUC__) \
 && __GNUC__ == 2 && __GNUC_MINOR__ == 96 && __GNUC_PATCHLEVEL__ == 0

#include <math.h>
namespace scitbx { namespace math {
  using ::atanh;
}} // namespace scitbx::math

#else

#include <boost/math/special_functions/atanh.hpp>
namespace scitbx { namespace math {
  using boost::math::atanh;
}} // namespace scitbx::math

#endif

#endif // SCITBX_MATH_ATANH_H
