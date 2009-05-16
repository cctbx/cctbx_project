#ifndef SCITBX_MATH_NUMERIC_LIMITS
#define SCITBX_MATH_NUMERIC_LIMITS

#include <limits>
#include <cmath>

namespace scitbx { namespace math {

/// Extension of std::numeric_limits for LAPACK-style floating point info
template <typename T>
class numeric_limits : public std::numeric_limits<T>
{
  public:
    /// Safe minimum, such that 1/safe_min does not overflow
    static T safe_min() {
      return std::pow(T(std::numeric_limits<T>::radix),
                      std::numeric_limits<T>::min_exponent);
    }

    /// Epsilon times radix
    static T epsilon_x_radix() {
      return std::numeric_limits<T>::epsilon()*std::numeric_limits<T>::radix;
    }
};


}} // scitbx::math

#endif // GUARD
