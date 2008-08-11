#ifndef SCITBX_MATH_ERF_ENGINE_H
#define SCITBX_MATH_ERF_ENGINE_H

#include <scitbx/math/erf/constants.h>
#include <cmath>

namespace scitbx { namespace math {

  //! Port of http://www.netlib.org/specfun/erf (as of 2003 Dec 03).
  /*! From the FORTRAN documentation:
<pre>
This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
for a real argument  x.  It contains three FUNCTION type
subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
and one SUBROUTINE type subprogram, CALERF.  The calling
statements for the primary entries are:

                Y=ERF(X)     (or   Y=DERF(X)),

                Y=ERFC(X)    (or   Y=DERFC(X)),
and
                Y=ERFCX(X)   (or   Y=DERFCX(X)).

The routine  CALERF  is intended for internal packet use only,
all computations within the packet being concentrated in this
routine.  The function subprograms invoke  CALERF  with the
statement

       CALL CALERF(ARG,RESULT,JINT)

where the parameter usage is as follows

   Function                     Parameters for CALERF
    call              ARG                  Result          JINT

  ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
  ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
  ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2

The main computation evaluates near-minimax approximations
from "Rational Chebyshev approximations for the error function"
by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
transportable program uses rational functions that theoretically
approximate  erf(x)  and  erfc(x)  to at least 18 significant
decimal digits.  The accuracy achieved depends on the arithmetic
system, the compiler, the intrinsic functions, and proper
selection of the machine-dependent constants.
</pre>
  */

  template <typename FloatType,
            typename AintIntType=long>
  struct erf_engine : erf_constants::machine_dependent<FloatType>
  {
    FloatType
    compute(FloatType const& arg, int jint)
    {
      FloatType sixten = 16.;
      FloatType four = 4.;
      FloatType one = 1.;
      FloatType half = .5;
      FloatType two = 2.;
      FloatType zero = 0.;
      FloatType sqrpi = 5.6418958354775628695e-1;
      FloatType thresh = 0.46875;
      FloatType a[] = {3.16112374387056560e00, 1.13864154151050156e02,
                       3.77485237685302021e02, 3.20937758913846947e03,
                       1.85777706184603153e-1
                      };
      FloatType b[] = {2.36012909523441209e01, 2.44024637934444173e02,
                       1.28261652607737228e03, 2.84423683343917062e03
                      };
      FloatType c[] = {5.64188496988670089e-1, 8.88314979438837594e0,
                       6.61191906371416295e01, 2.98635138197400131e02,
                       8.81952221241769090e02, 1.71204761263407058e03,
                       2.05107837782607147e03, 1.23033935479799725e03,
                       2.15311535474403846e-8
                      };
      FloatType d[] = {1.57449261107098347e01, 1.17693950891312499e02,
                       5.37181101862009858e02, 1.62138957456669019e03,
                       3.29079923573345963e03, 4.36261909014324716e03,
                       3.43936767414372164e03, 1.23033935480374942e03
                      };
      FloatType p[] = {3.05326634961232344e-1, 3.60344899949804439e-1,
                       1.25781726111229246e-1, 1.60837851487422766e-2,
                       6.58749161529837803e-4, 1.63153871373020978e-2
                      };
      FloatType q[] = {2.56852019228982242e00, 1.87295284992346047e00,
                       5.27905102951428412e-1, 6.05183413124413191e-2,
                       2.33520497626869185e-3
                      };
      FloatType result;
      FloatType x = arg;
      FloatType y = x;
      if (y < 0) y = -y;
      if (y <= thresh) {
            FloatType ysq = zero;
            if (y > this->xsmall) ysq = y * y;
            FloatType xnum = a[5-1]*ysq;
            FloatType xden = ysq;
            for(std::size_t i=1;i<=3;i++) {
               xnum = (xnum + a[i-1]) * ysq;
               xden = (xden + b[i-1]) * ysq;
            }
            result = x * (xnum + a[4-1]) / (xden + b[4-1]);
            if (jint != 0) result = one - result;
            if (jint == 2) result = exp_(ysq) * result;
            return result;
         }
         else if (y <= four) {
            FloatType xnum = c[9-1]*y;
            FloatType xden = y;
            for(std::size_t i=1;i<=7;i++) {
               xnum = (xnum + c[i-1]) * y;
               xden = (xden + d[i-1]) * y;
            }
            result = (xnum + c[8-1]) / (xden + d[8-1]);
            if (jint != 2) {
               FloatType ysq = aint_(y*sixten)/sixten;
               FloatType del_ = (y-ysq)*(y+ysq);
               result = exp_(-ysq*ysq) * exp_(-del_) * result;
            }
         }
         else {
            result = zero;
            if (y >= this->xbig) {
               if ((jint != 2) || (y >= this->xmax)) goto label_300;
               if (y >= this->xhuge) {
                  result = sqrpi / y;
                  goto label_300;
               }
            }
            FloatType ysq = one / (y * y);
            FloatType xnum = p[6-1]*ysq;
            FloatType xden = ysq;
            for(std::size_t i=1;i<=4;i++) {
               xnum = (xnum + p[i-1]) * ysq;
               xden = (xden + q[i-1]) * ysq;
            }
            result = ysq *(xnum + p[5-1]) / (xden + q[5-1]);
            result = (sqrpi -  result) / y;
            if (jint != 2) {
               FloatType ysq = aint_(y*sixten)/sixten;
               FloatType del_ = (y-ysq)*(y+ysq);
               result = exp_(-ysq*ysq) * exp_(-del_) * result;
            }
      }
      label_300:
      if (jint == 0) {
            result = (half - result) + half;
            if (x < zero) result = -result;
         }
         else if (jint == 1) {
            if (x < zero) result = two - result;
         }
         else {
            if (x < zero) {
               if (x < this->xneg) {
                     result = this->xinf;
                  }
                  else {
                     FloatType ysq = aint_(x*sixten)/sixten;
                     FloatType del_ = (x-ysq)*(x+ysq);
                     y = exp_(ysq*ysq) * exp_(del_);
                     result = (y+y) - result;
               }
            }
      }
      return result;
    }

    static FloatType
    exp_(FloatType const& x)
    {
      return static_cast<FloatType>(std::exp(x));
    }

    static FloatType
    aint_(FloatType const& x)
    {
      return static_cast<FloatType>(static_cast<AintIntType>(x));
    }
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_ERF_ENGINE_H
