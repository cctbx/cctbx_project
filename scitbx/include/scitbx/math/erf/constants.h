#ifndef SCITBX_MATH_ERF_CONSTANTS_H
#define SCITBX_MATH_ERF_CONSTANTS_H

#include <cstddef>

namespace scitbx { namespace math {

//! Port of http://www.netlib.org/specfun/erf (as of 2003 Dec 03).
namespace erf_constants {

  //! Port of http://www.netlib.org/specfun/erf (as of 2003 Dec 03).
  /*! From the FORTRAN documentation:
<pre>
Explanation of machine-dependent constants

  XMIN   = the smallest positive floating-point number.
  XINF   = the largest positive finite floating-point number.
  XNEG   = the largest negative argument acceptable to ERFCX;
           the negative of the solution to the equation
           2*exp(x*x) = XINF.
  XSMALL = argument below which erf(x) may be represented by
           2*x/sqrt(pi)  and above which  x*x  will not underflow.
           A conservative value is the largest machine number X
           such that   1.0 + X = 1.0   to machine precision.
  XBIG   = largest argument acceptable to ERFC;  solution to
           the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
           W(x) = exp(-x*x)/[x*sqrt(pi)].
  XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
           machine precision.  A conservative value is
           1/[2*sqrt(XSMALL)]
  XMAX   = largest acceptable argument to ERFCX; the minimum
           of XINF and 1/[sqrt(pi)*XMIN].
</pre>
   */
  template <typename FloatType>
  struct machine_dependent_base
  {
    machine_dependent_base(
      FloatType const& xinf_,
      FloatType const& xneg_,
      FloatType const& xsmall_,
      FloatType const& xbig_,
      FloatType const& xhuge_,
      FloatType const& xmax_)
    :
      xinf(xinf_),
      xneg(xneg_),
      xsmall(xsmall_),
      xbig(xbig_),
      xhuge(xhuge_),
      xmax(xmax_)
    {}

    FloatType xinf, xneg, xsmall, xbig, xhuge, xmax;
  };

  template <typename FloatType,
            std::size_t SizeOfFloatType=sizeof(FloatType)>
  struct machine_dependent;

  template <typename FloatType>
  struct machine_dependent<FloatType, 4> : machine_dependent_base<FloatType>
  {
    machine_dependent() : machine_dependent_base<FloatType>(
      3.40e38, -9.382, 5.96e-8, 9.194, 2.90e3, 4.79e37)
    {}
  };

  template <typename FloatType>
  struct machine_dependent<FloatType, 8> : machine_dependent_base<FloatType>
  {
    machine_dependent() : machine_dependent_base<FloatType>(
      1.79e308, -26.628e0, 1.11e-16, 26.543e0, 6.71e7, 2.53e307)
    {}
  };

}}} // namespace scitbx::math::erf_constants

#endif // SCITBX_MATH_ERF_CONSTANTS_H
