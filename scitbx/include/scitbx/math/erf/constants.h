#ifndef SCITBX_MATH_ERF_CONSTANTS_H
#define SCITBX_MATH_ERF_CONSTANTS_H

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
  struct machine_dependent
  {
    static const FloatType xinf;
    static const FloatType xneg;
    static const FloatType xsmall;
    static const FloatType xbig;
    static const FloatType xhuge;
    static const FloatType xmax;
  };

  template<> const float machine_dependent<float>::xinf = 3.40e+38;
  template<> const float machine_dependent<float>::xneg = -9.382e0;
  template<> const float machine_dependent<float>::xsmall = 5.96e-8;
  template<> const float machine_dependent<float>::xbig = 9.194e0;
  template<> const float machine_dependent<float>::xhuge = 2.90e3;
  template<> const float machine_dependent<float>::xmax = 4.79e37;

  template<> const double machine_dependent<double>::xinf = 1.79e308;
  template<> const double machine_dependent<double>::xneg = -26.628e0;
  template<> const double machine_dependent<double>::xsmall = 1.11e-16;
  template<> const double machine_dependent<double>::xbig = 26.543e0;
  template<> const double machine_dependent<double>::xhuge = 6.71e7;
  template<> const double machine_dependent<double>::xmax = 2.53e307;

}}} // namespace scitbx::math::erf_constants

#endif // SCITBX_MATH_ERF_CONSTANTS_H
