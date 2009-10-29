#ifndef SCITBX_MATH_APPROX_EQUAL_H
#define SCITBX_MATH_APPROX_EQUAL_H

#include <algorithm>
#include <cmath>
#include <scitbx/math/traits.h>
#include <limits>

namespace scitbx { namespace math {

  /** \brief A functor testing whether the relative error
    between two floating point values is smaller than some tolerance.

    This is recommended a method to compare two values which should only
    differ because of rounding errors. The generic type NumType may be
    a scalar floating point type, a complex number type or even some vector of
    floating point values, providing that std::abs(x) exists
    and that abs_traits<NumType> has the correct specialisation.

    The implementation first tests whether the absolute difference
    between the two values is very small, in which case, they are deemed to be
    equal to zero.
  */
  template <class NumType>
  struct approx_equal_relatively
  {
    typedef NumType num_type;
    typedef typename abs_traits<num_type>::result_type amplitude_type;

    amplitude_type relative_error, near_zero_threshold;

    approx_equal_relatively(amplitude_type relative_error) :
      relative_error(relative_error),
      near_zero_threshold(std::numeric_limits<amplitude_type>::min())
    {}

    bool operator()(num_type const &x, num_type const &y) const {
      amplitude_type d = std::abs(x - y);
      if (d <= near_zero_threshold) return true;
      amplitude_type a = std::abs(x), b = std::abs(y);
      if (d/std::max(a, b) <= relative_error) return true;
      return false;
    }
  };

  /** \brief A functor testing whether the relative error
    between two floating point values is smaller than some tolerance.

    The generic type NumType may be
    a scalar floating point type, a complex number type or even some vector of
    floating point values, providing that std::abs(x) exists
    and that abs_traits<NumType> has the correct specialisation.
   */
  template <class NumType>
  struct approx_equal_absolutely
  {
    typedef NumType num_type;
    typedef typename abs_traits<num_type>::result_type amplitude_type;

    amplitude_type absolute_error;

    approx_equal_absolutely(amplitude_type absolute_error)
      : absolute_error(absolute_error)
    {}

    bool operator()(num_type const &x, num_type const &y) const {
      return std::abs(x - y) <= absolute_error;
    }
  };
}}

#endif
