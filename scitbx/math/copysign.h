#ifndef SCITBX_MATH_COPY_SIGN_H
#define SCITBX_MATH_COPY_SIGN_H

namespace scitbx { namespace math {

  /** \fn double copysign(double x, double y)
       \brief  Return x with its sign changed to match the sign of y.

       \fn float copysign(float x, float y)
       \copydoc double copysign(double x, double y)
  */

}}

#if defined(__GNUC__) && __GNUC__ == 4
  #include <math.h>
  namespace scitbx { namespace math {
    inline
    double copysign(double x, double y) {
      return ::copysign(x,y);
    }
    inline
    float copysign(float x, float y) {
      return ::copysignf(x,y);
    }
  }}
#elif defined(BOOST_MSVC) && (BOOST_MSVC >= 1310)
  // VC++ 7.1 or newer
  #include <math.h>
  #include <boost/math/special_functions/sign.hpp>
  namespace scitbx { namespace math {
    // VC++ doesn't provide the float version
    inline
    double copysign(double x, double y) {
      return ::_copysign(x,y);
    }
    inline
    float copysign(float x, float y) {
      return boost::math::copysign(x,y);
    }
  }}
#else
  #include <boost/math/special_functions/sign.hpp>
  namespace scitbx { namespace math {
    using boost::math::copysign;
  }}
#endif

#endif // GUARD
