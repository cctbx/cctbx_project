#ifndef SCITBX_MATH_ERF_H
#define SCITBX_MATH_ERF_H

#include <scitbx/math/erf/engine.h>

namespace scitbx { namespace math {

  //! Approximate values for erf(x).
  /*! See also: scitbx::math::erf_engine
   */
  template <typename FloatType>
  inline
  FloatType
  erf(FloatType const& x)
  {
    return erf_engine<FloatType>().compute(x, 0);
  }

  //! Approximate values for erfc(x).
  /*! See also: scitbx::math::erf_engine
   */
  template <typename FloatType>
  inline
  FloatType
  erfc(FloatType const& x)
  {
    return erf_engine<FloatType>().compute(x, 1);
  }

  //! Approximate values for exp(x*x) * erfc(x).
  /*! See also: scitbx::math::erf_engine
   */
  template <typename FloatType>
  inline
  FloatType
  erfcx(FloatType const& x)
  {
    return erf_engine<FloatType>().compute(x, 2);
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_ERF_H
