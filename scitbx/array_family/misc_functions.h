#ifndef SCITBX_ARRAY_FAMILY_MISC_FUNCTIONS_H
#define SCITBX_ARRAY_FAMILY_MISC_FUNCTIONS_H

#include <boost/type_traits/is_unsigned.hpp>
#include <cmath>
#include <cstdlib>

namespace scitbx { namespace fn {

  namespace details {

    template <typename NumType>
    class absolute
    {
    public:
      NumType operator()(NumType x) {
        return (*this)(x, typename boost::is_unsigned<NumType>::type());
      }

    private:
      /* boost::is_unsigned<my_type> is false by default, which makes
         this version the default, and the other one the version used only
         for those types explicitly marked as unsigned.
       */
      NumType operator()(NumType x, boost::false_type is_unsigned) {
        return std::abs(x);
      }

      NumType operator()(NumType x, boost::true_type is_unsigned) {
        return x;
      }
    };
  }

  //! Absolute value.
  template <typename NumType>
  inline
  NumType
  absolute(NumType const& x)
  {
    return details::absolute<NumType>()(x);
  }

  //! Floating-point modulus with strictly positive result.
  template <typename FloatType1, typename FloatType2>
  inline
  FloatType1
  fmod_positive(FloatType1 const& x, FloatType2 const& y)
  {
    FloatType1 result = std::fmod(x, y);
    while (result < FloatType1(0)) result += y;
    return result;
  }

  //! Square of x.
  template <typename NumType>
  inline
  NumType
  pow2(NumType const& x)
  {
    return x * x;
  }

  //! Cube of x.
  // Not available as array function.
  template <typename NumType>
  inline
  NumType
  pow3(NumType const& x)
  {
    return x * x * x;
  }

  //! Fourth power of x.
  // Not available as array function.
  template <typename NumType>
  inline
  NumType
  pow4(NumType const& x)
  {
    return pow2(pow2(x));
  }

  //! Tests if abs(a-b) <= tolerance.
  template <class FloatType>
  bool
  approx_equal(FloatType const& a,
               FloatType const& b,
               FloatType const& tolerance)
  {
    FloatType diff = a - b;
    if (diff < 0.) diff = -diff;
    if (diff <= tolerance) return true;
    return false;
  }

  //! Helper function object for array operations.
  template <typename ResultType,
            typename ArgumentType>
  struct functor_absolute
  {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType const& x) const {
      return ResultType(absolute(x));
    }
  };

  //! Helper function object for array operations.
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2>
  struct functor_fmod_positive
  {
    typedef ResultType result_type;
    ResultType
    operator()(ArgumentType1 const& x, ArgumentType2 const& y) const
    {
      return ResultType(fmod_positive(x, y));
    }
  };

  //! Helper function object for array operations.
  template <typename ResultType,
            typename ArgumentType>
  struct functor_pow2
  {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType const& x) const {
      return ResultType(pow2(x));
    }
  };

  //! Helper function object for array operations.
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2,
            typename ArgumentType3>
  struct functor_approx_equal
  {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType1 const& x,
                          ArgumentType2 const& y,
                          ArgumentType3 const& z) const {
    return ResultType(approx_equal(x, y, z)); }
  };

}} // namespace scitbx::fn

#endif // SCITBX_ARRAY_FAMILY_MISC_FUNCTIONS_H
