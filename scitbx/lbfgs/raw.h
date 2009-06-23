// The help of Objexx in converting this file is greatly appreciated!
// http://objexx.com/

#ifndef SCITBX_LBFGS_RAW_H
#define SCITBX_LBFGS_RAW_H

#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cmath>

namespace scitbx { namespace lbfgs { namespace raw {

  //! Emulation of 1-dimensional FORTRAN arrays with offset 1.
  template <typename ElementType>
  class const_ref1
  {
    protected:
      const ElementType* begin_;
      int n_;

    public:
      const_ref1() {}

      const_ref1(const ElementType* begin, int n) : begin_(begin), n_(n) {}

      const ElementType*
      begin() const { return begin_; }

      const ElementType*
      end() const { return begin_+n_; }

      const ElementType&
      operator()(int i) const { return begin_[i-1]; }

      const_ref1
      get1(int i, int n) const { return const_ref1(begin_+i-1, n); }
  };

  template <typename ElementType>
  class ref1 : public const_ref1<ElementType>
  {
    public:
      ref1() {}

      ref1(ElementType* begin, int n) : const_ref1<ElementType>(begin, n) {}

      ElementType*
      begin() const { return const_cast<ElementType*>(this->begin_); }

      ElementType*
      end() const { return this->begin()+this->n_; }

      ElementType&
      operator()(int i) const { return this->begin()[i-1]; }

      ref1
      get1(int i, int n) const { return ref1(this->begin()+i-1, n); }
  };

  inline
  void
  daxpy(
    int const n,
    double const da,
    const_ref1<double> const& dx,
    int const incx,
    ref1<double> const& dy,
    int const incy)
  {
    // constant times a vector plus a vector.
    // uses unrolled loops for increments equal to one.
    // jack dongarra, linpack, 3/11/78.f

    int ix, iy, m, mp1;

    if ( n <= 0 ) return;
    if ( da == 0.0 ) return;
    if ( incx == 1 && incy == 1 ) goto L20;

    // code for unequal increments or equal increments
    // not equal to 1

    ix = 1;
    iy = 1;
    if ( incx < 0 ) ix = ( -n + 1)*incx + 1;
    if ( incy < 0 ) iy = ( -n + 1)*incy + 1;
    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      dy(iy) += da*dx(ix);
      ix += incx;
      iy += incy;
    }
    return;

    // code for both increments equal to 1
    //
    // clean-up loop

  L20:
    m = n % 4;
    if ( m == 0 ) goto L40;
    for ( int i = 1, i_end = m; i <= i_end; ++i ) {
      dy(i) += da*dx(i);
    }
    if ( n < 4 ) return;
  L40:
    mp1 = m + 1;
    for ( int i = mp1, i_end = n; i <= i_end; i += 4 ) {
      dy(i) += da*dx(i);
      dy(i + 1) += da*dx(i + 1);
      dy(i + 2) += da*dx(i + 2);
      dy(i + 3) += da*dx(i + 3);
    }
  }

  inline
  double
  ddot(
    int const n,
    const_ref1<double> const& dx,
    int const incx,
    const_ref1<double> const& dy,
    int const incy)
  {
    double ddot; // Return value

    // forms the dot product of two vectors.
    // uses unrolled loops for increments equal to one.
    // jack dongarra, linpack, 3/11/78.f

    double dtemp;
    int ix, iy, m, mp1;

    ddot = 0.0;
    dtemp = 0.0;
    if ( n <= 0 ) return  ddot;
    if ( incx == 1 && incy == 1 ) goto L20;

    // code for unequal increments or equal increments
    // not equal to 1

    ix = 1;
    iy = 1;
    if ( incx < 0 ) ix = ( -n + 1)*incx + 1;
    if ( incy < 0 ) iy = ( -n + 1)*incy + 1;
    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      dtemp += dx(ix)*dy(iy);
      ix += incx;
      iy += incy;
    }
    ddot = dtemp;
    return  ddot;

    //        code for both increments equal to 1
    //
    //        clean-up loop

  L20:
    m = n % 5;
    if ( m == 0 ) goto L40;
    for ( int i = 1, i_end = m; i <= i_end; ++i ) {
      dtemp += dx(i)*dy(i);
    }
    if ( n < 5 ) goto L60;
  L40:
    mp1 = m + 1;
    for ( int i = mp1, i_end = n; i <= i_end; i += 5 ) {
      dtemp += dx(i)*dy(i)
             + dx(i + 1)*dy(i + 1)
             + dx(i + 2)*dy(i + 2)
             + dx(i + 3)*dy(i + 3)
             + dx(i + 4)*dy(i + 4);
    }
  L60:
    ddot = dtemp;
    return dtemp;
  }

  inline
  void
  mcstep(
    double & stx,
    double & fx,
    double & dx,
    double & sty,
    double & fy,
    double & dy,
    double & stp,
    double const fp,
    double const dp,
    bool & brackt,
    double const stpmin,
    double const stpmax,
    int & info)
  {
    bool bound;

    //     SUBROUTINE MCSTEP
    //
    //     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
    //     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
    //     A MINIMIZER OF THE FUNCTION.
    //
    //     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
    //     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
    //     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
    //     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
    //     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
    //     WITH ENDPOINTS STX AND STY.
    //
    //     THE SUBROUTINE STATEMENT IS
    //
    //       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
    //                        STPMIN,STPMAX,INFO)
    //
    //     WHERE
    //
    //       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
    //         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
    //         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
    //         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
    //         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
    //
    //       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
    //         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
    //         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
    //         UPDATED APPROPRIATELY.
    //
    //       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
    //         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
    //         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
    //         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
    //
    //       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
    //         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
    //         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
    //         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
    //
    //       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
    //         AND UPPER BOUNDS FOR THE STEP.
    //
    //       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    //         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
    //         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
    //         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
    //
    //     SUBPROGRAMS CALLED
    //
    //       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
    //
    //     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
    //     JORGE J. MORE', DAVID J. THUENTE
    //
    double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
    info = 0;

    //     CHECK THE INPUT PARAMETERS FOR ERRORS.

    if ( (brackt && (stp <= std::min(stx,sty) || stp >= std::max(stx, sty)))
         || dx*(stp - stx) >= 0.0
         || stpmax < stpmin ) return;

    //     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.

    sgnd = dp*(dx/std::abs(dx));

    //     FIRST CASE. A HIGHER FUNCTION VALUE.
    //     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
    //     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
    //     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.

    if ( fp > fx ) {
      info = 1;
      bound = true;
      theta = 3*(fx - fp)/(stp - stx) + dx + dp;
      s = std::max(std::abs(theta),std::max(std::abs(dx),std::abs(dp)));
      gamma = s*std::sqrt(std::pow( (theta/s), 2 ) - (dx/s)*(dp/s));
      if ( stp < stx ) gamma = -gamma;
      p = (gamma - dx) + theta;
      q = ((gamma - dx) + gamma) + dp;
      r = p/q;
      stpc = stx + r*(stp - stx);
      stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/2)*(stp - stx);
      if ( std::abs(stpc - stx) < std::abs(stpq - stx) ) {
        stpf = stpc;
      }
      else {
        stpf = stpc + (stpq - stpc)/2;
      }
      brackt = true;

    //     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
    //     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
    //     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
    //     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.

    }
    else if ( sgnd < 0.0 ) {
      info = 2;
      bound = false;
      theta = 3*(fx - fp)/(stp - stx) + dx + dp;
      s = std::max(std::abs(theta),std::max(std::abs(dx),std::abs(dp)));
      gamma = s*std::sqrt(std::pow( (theta/s), 2 ) - (dx/s)*(dp/s));
      if ( stp > stx ) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dx;
      r = p/q;
      stpc = stp + r*(stx - stp);
      stpq = stp + (dp/(dp - dx))*(stx - stp);
      if ( std::abs(stpc - stp) > std::abs(stpq - stp) ) {
        stpf = stpc;
      }
      else {
        stpf = stpq;
      }
      brackt = true;

    //     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
    //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
    //     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
    //     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
    //     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
    //     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
    //     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
    //     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.

    }
    else if ( std::abs(dp) < std::abs(dx) ) {
      info = 3;
      bound = true;
      theta = 3*(fx - fp)/(stp - stx) + dx + dp;
      s = std::max(std::abs(theta),std::max(std::abs(dx),std::abs(dp)));

    //        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
    //        TO INFINITY IN THE DIRECTION OF THE STEP.

      gamma = s*std::sqrt(std::max(0.0,
        std::pow( (theta/s), 2 ) - (dx/s)*(dp/s)));
      if ( stp > stx ) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = (gamma + (dx - dp)) + gamma;
      r = p/q;
      if ( r < 0.0 && gamma != 0.0 ) {
        stpc = stp + r*(stx - stp);
      }
      else if ( stp > stx ) {
        stpc = stpmax;
      }
      else {
        stpc = stpmin;
      }
      stpq = stp + (dp/(dp - dx))*(stx - stp);
      if ( brackt ) {
        if ( std::abs(stp - stpc) < std::abs(stp - stpq) ) {
          stpf = stpc;
        }
        else {
          stpf = stpq;
        }
      }
      else {
        if ( std::abs(stp - stpc) > std::abs(stp - stpq) ) {
          stpf = stpc;
        }
        else {
          stpf = stpq;
        }
      }

    //     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
    //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
    //     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
    //     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.

    }
    else {
      info = 4;
      bound = false;
      if ( brackt ) {
        theta = 3*(fp - fy)/(sty - stp) + dy + dp;
        s = std::max(std::abs(theta),std::max(std::abs(dy),std::abs(dp)));
        gamma = s*std::sqrt(std::pow( (theta/s), 2 ) - (dy/s)*(dp/s));
        if ( stp > sty ) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dy;
        r = p/q;
        stpc = stp + r*(sty - stp);
        stpf = stpc;
      }
      else if ( stp > stx ) {
        stpf = stpmax;
      }
      else {
        stpf = stpmin;
      }
    }

    //     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
    //     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.

    if ( fp > fx ) {
      sty = stp;
      fy = fp;
      dy = dp;
    }
    else {
      if ( sgnd < 0.0 ) {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }

    //     COMPUTE THE NEW STEP AND SAFEGUARD IT.

    stpf = std::min(stpmax,stpf);
    stpf = std::max(stpmin,stpf);
    stp = stpf;
    if ( brackt && bound ) {
      if ( sty > stx ) {
        stp = std::min(stx + 0.66*(sty - stx),stp);
      }
      else {
        stp = std::max(stx + 0.66*(sty - stx),stp);
      }
    }
  }

  namespace LB3 {
    static int lp = 6;
    static double gtol = 9.0e-01;
    static double stpmin = 1.0e-20;
    static double stpmax = 1.0e+20;
  }

  inline
  void
  mcsrch(
    int const n,
    ref1<double> const& x,
    double const f,
    const_ref1<double> const& g,
    const_ref1<double> const& s,
    double & stp,
    double const ftol,
    double const xtol,
    int const maxfev,
    int & info,
    int & nfev,
    ref1<double> const& wa)
  {
    using namespace LB3;

    //                     SUBROUTINE MCSRCH
    //
    //     A slight modification of the subroutine CSRCH of More' and Thuente.
    //     The changes are to allow reverse communication, and do not affect
    //     the performance of the routine.
    //
    //     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
    //     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
    //
    //     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
    //     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
    //     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
    //     MINIMIZER OF THE MODIFIED FUNCTION
    //
    //          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
    //
    //     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
    //     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
    //     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
    //     CONTAINS A MINIMIZER OF F(X+STP*S).
    //
    //     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
    //     THE SUFFICIENT DECREASE CONDITION
    //
    //           F(X+STP*S) <= F(X) + FTOL*STP*(GRADF(X)'S),
    //
    //     AND THE CURVATURE CONDITION
    //
    //           std::abs(GRADF(X+STP*S)'S)) <= GTOL*std::abs(GRADF(X)'S).
    //
    //     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
    //     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
    //     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
    //     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
    //     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
    //     SATISFIES THE SUFFICIENT DECREASE CONDITION.
    //
    //     THE SUBROUTINE STATEMENT IS
    //
    //        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
    //     WHERE
    //
    //       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
    //         OF VARIABLES.
    //
    //       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
    //         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
    //         X + STP*S.
    //
    //       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
    //         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
    //
    //       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
    //         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
    //         OF F AT X + STP*S.
    //
    //       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
    //         SEARCH DIRECTION.
    //
    //       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
    //         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
    //         STP CONTAINS THE FINAL ESTIMATE.
    //
    //       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
    //         communication implementation GTOL is defined in a COMMON
    //         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
    //         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
    //         SATISFIED.
    //
    //       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
    //         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
    //         IS AT MOST XTOL.
    //
    //       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
    //         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
    //         communication implementatin they are defined in a COMMON
    //         statement).
    //
    //       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
    //         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
    //         MAXFEV BY THE END OF AN ITERATION.
    //
    //       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    //
    //         INFO = 0  IMPROPER INPUT PARAMETERS.
    //
    //         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
    //
    //         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
    //                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
    //
    //         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
    //                   IS AT MOST XTOL.
    //
    //         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
    //
    //         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
    //
    //         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
    //
    //         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
    //                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
    //                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
    //                   TOLERANCES MAY BE TOO SMALL.
    //
    //       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
    //         CALLS TO FCN.
    //
    //       WA IS A WORK ARRAY OF LENGTH N.
    //
    //     SUBPROGRAMS CALLED
    //
    //       MCSTEP
    //
    //       FORTRAN-SUPPLIED...ABS,MAX,MIN
    //
    //     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
    //     JORGE J. MORE', DAVID J. THUENTE
    //
    //     **********
    static int infoc;
    static bool brackt, stage1;
    static double
      dg, dgm, dginit, dgtest, dgx, dgxm, dgy, dgym, finit, ftest1,
      fm, fx, fxm, fy, fym, stx, sty, stmin, stmax, width, width1;
    static double p5 = 0.5;
    static double p66 = 0.66;
    static double xtrapf = 4.0;
    static double zero = 0.0;

    if ( info == -1 ) goto L45;
    infoc = 1;

    //     CHECK THE INPUT PARAMETERS FOR ERRORS.

    if ( n <= 0 || stp <= zero || ftol < zero || gtol < zero || xtol < zero
         || stpmin < zero || stpmax < stpmin || maxfev <= 0 ) return;

    //     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
    //     AND CHECK THAT S IS A DESCENT DIRECTION.

    dginit = zero;
    for ( int j = 1, j_end = n; j <= j_end; ++j ) {
      dginit += g(j) * s(j);
    }
    if ( dginit >= zero ) {
      std::printf("\n  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION\n");
      return;
    }

    //     INITIALIZE LOCAL VARIABLES.

    brackt = false;
    stage1 = true;
    nfev = 0;
    finit = f;
    dgtest = ftol*dginit;
    width = stpmax - stpmin;
    width1 = width/p5;
    std::copy(x.begin(), x.end(), wa.begin());

    //     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
    //     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
    //     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
    //     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
    //     THE INTERVAL OF UNCERTAINTY.
    //     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
    //     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.

    stx = zero;
    fx = finit;
    dgx = dginit;
    sty = zero;
    fy = finit;
    dgy = dginit;

    //     START OF ITERATION.

  L30:

    //     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
    //     TO THE PRESENT INTERVAL OF UNCERTAINTY.

    if ( brackt ) {
      stmin = std::min(stx,sty);
      stmax = std::max(stx,sty);
    }
    else {
      stmin = stx;
      stmax = stp + xtrapf*(stp - stx);
    }

    //     FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.

    stp = std::max(stp,stpmin);
    stp = std::min(stp,stpmax);

    //     IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
    //     STP BE THE LOWEST POINT OBTAINED SO FAR.

    if ( (brackt && (stp <= stmin || stp >= stmax)) || nfev >= maxfev - 1 || infoc == 0 || (brackt && stmax - stmin <= xtol*stmax) ) stp = stx;

    //     EVALUATE THE FUNCTION AND GRADIENT AT STP
    //     AND COMPUTE THE DIRECTIONAL DERIVATIVE.
    //     We return to main program to obtain F and G.

    for ( int j = 1, j_end = n; j <= j_end; ++j ) {
      x(j) = wa(j) + stp*s(j);
    }
    info = -1;
    return;

  L45:
    info = 0;
    ++nfev;
    dg = zero;
    for ( int j = 1, j_end = n; j <= j_end; ++j ) {
      dg += g(j)*s(j);
    }
    ftest1 = finit + stp*dgtest;

    //     TEST FOR CONVERGENCE.

    if ( (brackt && (stp <= stmin || stp >= stmax)) || infoc == 0 ) info = 6;
    if ( stp == stpmax && f <= ftest1 && dg <= dgtest ) info = 5;
    if ( stp == stpmin && (f > ftest1 || dg >= dgtest) ) info = 4;
    if ( nfev >= maxfev ) info = 3;
    if ( brackt && stmax - stmin <= xtol*stmax ) info = 2;
    if ( f <= ftest1 && std::abs(dg) <= gtol*( - dginit) ) info = 1;

    //     CHECK FOR TERMINATION.

    if ( info != 0 ) return;

    //     IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
    //     FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.

    if ( stage1 && f <= ftest1 && dg >= std::min(ftol,gtol) *dginit ) stage1 = false;

    //     A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
    //     WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
    //     FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
    //     DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
    //     OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.

    if ( stage1 && f <= fx && f > ftest1 ) {

    //        DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.

      fm = f - stp*dgtest;
      fxm = fx - stx*dgtest;
      fym = fy - sty*dgtest;
      dgm = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

    //        CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
    //        AND TO COMPUTE THE NEW STEP.

      mcstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,stmin,stmax,infoc);

    //        RESET THE FUNCTION AND GRADIENT VALUES FOR F.

      fx = fxm + stx*dgtest;
      fy = fym + sty*dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    }
    else {

    //        CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
    //        AND TO COMPUTE THE NEW STEP.

      mcstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,stmin,stmax,infoc);
    }

    //     FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
    //     INTERVAL OF UNCERTAINTY.

    if ( brackt ) {
      if ( std::abs(sty - stx) >= p66*width1 ) stp = stx + p5*(sty - stx);
      width1 = width;
      width = std::abs(sty - stx);
    }

    //     END OF ITERATION.

    goto L30;
  }

  inline
  void
  lb1_show_vector(
    const_ref1<double> const& v,
    int const n)
  {
    for(int i=1;i<=n;i++) {
      std::printf("  %10.3E", v(i));
      if (i % 6 == 0 || i == n) std::printf("\n");
    }
  }

  inline
  void
  lb1(
    const_ref1<int> const& iprint,
    int const iter,
    int const nfun,
    double const gnorm,
    int const n,
    int const m,
    const_ref1<double> const& x,
    double const f,
    const_ref1<double> const& g,
    double const stp,
    bool const finish)
  {
    //     -------------------------------------------------------------
    //     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
    //     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
    //     -------------------------------------------------------------

    using namespace LB3;

    if ( iter == 0 ) {
      std::printf("*************************************************\n");
      std::printf(
        "  N=%5d   NUMBER OF CORRECTIONS=%2d\n       INITIAL VALUES\n", n, m);
      std::printf(
        " F= %10.3E   GNORM= %10.3E\n", f, gnorm);
      if ( iprint(2) >= 1 ) {
        std::printf(" VECTOR X= \n");
        lb1_show_vector(x, n);
        std::printf(" GRADIENT VECTOR G= \n");
        lb1_show_vector(g, n);
      }
      std::printf("*************************************************\n");
      std::printf("\n   I   NFN    FUNC        GNORM       STEPLENGTH\n\n");
    }
    else {
      if ( (iprint(1) == 0) && (iter != 1 && ! finish) ) return;
      if ( iprint(1) != 0 ) {
        if ( (iter - 1) % iprint(1) == 0 || finish ) {
          if ( iprint(2) > 1 && iter > 1 ) {
            std::printf(
              "\n   I   NFN    FUNC        GNORM       STEPLENGTH\n\n");
          }
          std::printf(
            "%4d %4d    %10.3E  %10.3E  %10.3E\n", iter, nfun, f, gnorm, stp);
        }
        else {
          return;
        }
      }
      else {
        if ( iprint(2) > 1 && finish ) {
          std::printf(
            "\n   I   NFN    FUNC        GNORM       STEPLENGTH\n\n");
        }
        std::printf(
          "%4d %4d    %10.3E  %10.3E  %10.3E\n", iter, nfun, f, gnorm, stp);
      }
      if ( iprint(2) == 2 || iprint(2) == 3 ) {
        if ( finish ) {
          std::printf(" FINAL POINT X= \n");
        }
        else {
          std::printf(" VECTOR X= \n");
        }
        lb1_show_vector(x, n);
        if ( iprint(2) == 3 ) {
          std::printf(" GRADIENT VECTOR G= \n");
          lb1_show_vector(g, n);
        }
      }
      if ( finish ) {
        std::printf(
          "\n THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.\n"
          " IFLAG = 0\n");
      }
    }

  }

struct lbfgs {

  double gnorm, stp1, ftol, stp, ys, sq, yr, beta, xnorm;
  int
    iter, nfun, point, ispt, iypt, maxfev, info, bound, npt, cp, nfev,
    inmc, iycn, iscn;
  bool finish;
  int prev_iflag;
  ref1<double> current_search_direction;

  lbfgs()
  :
    gnorm(-1e20), stp1(-1e20), ftol(-1e20), stp(-1e20), ys(-1e20), sq(-1e20),
    yr(-1e20), beta(-1e20), xnorm(-1e20),
    iter(-9999), nfun(-9999), point(-9999), ispt(-9999), iypt(-9999),
    maxfev(-9999), info(-9999), bound(-9999), npt(-9999), cp(-9999),
    nfev(-9999), inmc(-9999), iycn(-9999), iscn(-9999),
    finish(true),
    prev_iflag(-9999),
    current_search_direction(0,0)
  {}

  void
  operator()(
    int const n,
    int const m,
    ref1<double> const& x,
    double const f,
    const_ref1<double> const& g,
    int const diagco,
    ref1<double> const& diag,
    const_ref1<int> const& iprint,
    double const eps,
    double const xtol,
    ref1<double> const& w,
    int & iflag)
  {
    //        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
    //                          JORGE NOCEDAL
    //                        *** July 1990 ***
    //
    //
    //     This subroutine solves the unconstrained minimization problem
    //
    //                      min F(x),    x= (x1,x2,...,xN),
    //
    //      using the limited memory BFGS method. The routine is especially
    //      effective on problems involving a large number of variables. In
    //      a typical iteration of this method an approximation Hk to the
    //      inverse of the Hessian is obtained by applying M BFGS updates to
    //      a diagonal matrix Hk0, using information from the previous M steps.
    //      The user specifies the number M, which determines the amount of
    //      storage required by the routine. The user may also provide the
    //      diagonal matrices Hk0 if not satisfied with the default choice.
    //      The algorithm is described in "On the limited memory BFGS method
    //      for large scale optimization", by D. Liu and J. Nocedal,
    //      Mathematical Programming B 45 (1989) 503-528.f
    //
    //      The user is required to calculate the function value F and its
    //      gradient G. In order to allow the user complete control over
    //      these computations, reverse  communication is used. The routine
    //      must be called repeatedly under the control of the parameter
    //      IFLAG.
    //
    //      The steplength is determined at each iteration by means of the
    //      line search routine MCVSRCH, which is a slight modification of
    //      the routine CSRCH written by More' and Thuente.
    //
    //      The calling statement is
    //
    //          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
    //
    //      where
    //
    //     N       is an INTEGER variable that must be set by the user to the
    //             number of variab^les. It is not altered by the routine.
    //             Restriction: N>0.
    //
    //     M       is an INTEGER variable that must be set by the user to
    //             the number of corrections used in the BFGS update. It
    //             is not altered by the routine. Values of M less than 3 are
    //             not recommended; large values of M will result in excessive
    //             computing time. 3<= M <=7 is recommended. Restriction: M>0.
    //
    //     X       is a DOUBLEPRECISION array of length N. On initial entry
    //             it must be set by the user to the values of the initial
    //             estimate of the solution vector. On exit with IFLAG=0, it
    //             contains the values of the variables at the best point
    //             found (usually a solution).
    //
    //     F       is a DOUBLEPRECISION variable. Before initial entry and on
    //             a re-entry with IFLAG=1, it must be set by the user to
    //             contain the value of the function F at the point X.
    //
    //     G       is a DOUBLEPRECISION array of length N. Before initial
    //             entry and on a re-entry with IFLAG=1, it must be set by
    //             the user to contain the components of the gradient G at
    //             the point X.
    //
    //     DIAGCO  is a LOGICAL variable that must be set to true if the
    //             user  wishes to provide the diagonal matrix Hk0 at each
    //             iteration. Otherwise it should be set to false, in which
    //             case  LBFGS will use a default value described below. If
    //             DIAGCO is set to true the routine will return at each
    //             iteration of the algorithm with IFLAG=2, and the diagonal
    //              matrix Hk0  must be provided in the array DIAG.
    //
    //
    //     DIAG    is a DOUBLEPRECISION array of length N. If DIAGCO=true,
    //             then on initial entry or on re-entry with IFLAG=2, DIAG
    //             it must be set by the user to contain the values of the
    //             diagonal matrix Hk0.  Restriction: all elements of DIAG
    //             must be positive.
    //
    //     IPRINT  is an INTEGER array of length two which must be set by the
    //             user.
    //
    //             IPRINT(1) specifies the frequency of the output:
    //                IPRINT(1) < 0 : no output is generated,
    //                IPRINT(1) = 0 : output only at first and last iteration,
    //                IPRINT(1) > 0 : output every IPRINT(1) iterations.
    //
    //             IPRINT(2) specifies the type of output generated:
    //                IPRINT(2) = 0 : iteration count, number of function
    //                                evaluations, function value, norm of the
    //                                gradient, and steplength,
    //                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
    //                                variables and  gradient vector at the
    //                                initial point,
    //                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
    //                                variables,
    //                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
    //
    //
    //     EPS     is a positive DOUBLEPRECISION variable that must be set by
    //             the user, and determines the accuracy with which the solution
    //             is to be found. The subroutine terminates when
    //
    //                         ||G|| < EPS max(1,||X||),
    //
    //             where ||.|| denotes the Euclidean norm.
    //
    //     XTOL    is a  positive DOUBLEPRECISION variable that must be set by
    //             the user to an estimate of the machine precision (e.g.
    //             std::pow( 10, (-16) ) on a SUN station 3/60). The line
    //             search routine will terminate if the relative width of
    //             the interval of uncertainty is less than XTOL.
    //
    //     W       is a DOUBLEPRECISION array of length N(2M+1)+2M used as
    //             workspace for LBFGS. This array must not be altered by the
    //             user.
    //
    //     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
    //             to the subroutine. A return with IFLAG<0 indicates an error,
    //             and IFLAG=0 indicates that the routine has terminated without
    //             detecting errors. On a return with IFLAG=1, the user must
    //             evaluate the function F and gradient G. On a return with
    //             IFLAG=2, the user must provide the diagonal matrix Hk0.
    //
    //             The following negative values of IFLAG, detecting an error,
    //             are possible:
    //
    //              IFLAG=-1  The line search routine MCSRCH failed. The
    //                        parameter INFO provides more detailed information
    //                        (see also the documentation of MCSRCH):
    //
    //                       INFO = 0  IMPROPER INPUT PARAMETERS.
    //
    //                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
    //                                 UNCERTAINTY IS AT MOST XTOL.
    //
    //                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
    //                                 REQUIRED AT THE PRESENT ITERATION.
    //
    //                       INFO = 4  THE STEP IS TOO SMALL.
    //
    //                       INFO = 5  THE STEP IS TOO LARGE.
    //
    //                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
    //                                 THERE MAY NOT BE A STEP WHICH SATISFIES
    //                                 THE SUFFICIENT DECREASE AND CURVATURE
    //                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
    //
    //
    //              IFLAG=-2  The i-th diagonal element of the diagonal inverse
    //                        Hessian approximation, given in DIAG, is not
    //                        positive.
    //
    //              IFLAG=-3  Improper input parameters for LBFGS (N or M are
    //                        not positive).
    //
    //
    //
    //    ON THE DRIVER:
    //
    //    The program that calls LBFGS must contain the declaration:
    //
    //                       EXTERNAL LB2
    //
    //    LB2 is a BLOCK DATA that defines the default values of several
    //    parameters described in the COMMON section.
    //
    //
    //
    //    COMMON:
    //
    //     The subroutine contains one common area, which the user may wish to
    //    reference:
    //
    using namespace LB3;
    //
    //    MP  is an INTEGER variable with default value 6.f It is used as the
    //        unit number for the printing of the monitoring information
    //        controlled by IPRINT.
    //
    //    LP  is an INTEGER variable with default value 6.f It is used as the
    //        unit number for the printing of error messages. This printing
    //        may be suppressed by setting LP to a non-positive value.
    //
    //    GTOL is a DOUBLEPRECISION variable with default value 0.9, which
    //        controls the accuracy of the line search routine MCSRCH. If the
    //        function and gradient evaluations are inexpensive with respect
    //        to the cost of the iteration (which is sometimes the case when
    //        solving very large problems) it may be advantageous to set GTOL
    //        to a small value. A typical small value is 0.1.  Restriction:
    //        GTOL should be greater than 1.e-04.f
    //
    //    STPMIN and STPMAX are non-negative DOUBLEPRECISION variables which
    //        specify lower and uper bounds for the step in the line search.
    //        Their default values are 1.e-20 and 1.e+20, respectively. These
    //        values need not be modified unless the exponents are too large
    //        for the machine being used, or unless the problem is extremely
    //        badly scaled (in which case the exponents should be increased).
    //
    //
    //  MACHINE DEPENDENCIES
    //
    //        The only variables that are machine-dependent are XTOL,
    //        STPMIN and STPMAX.
    //
    //
    //  GENERAL INFORMATION
    //
    //    Other routines called directly:  DAXPY, ddot, lb1, mcsrch
    //
    //    Input/Output  :  No input; diagnostic messages on unit MP and
    //                     error messages on unit LP.
    //
    //
    //     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    const static double one = 1.0;
    const static double zero = 0.0;

    if (diagco < 0 || diagco > 2) {
      throw std::runtime_error("raw::lbfgs: invalid diagco value.");
    }

    static const char* iflag_minus_2_format =
      "\n IFLAG= -2\n THE%5d-TH DIAGONAL ELEMENT OF THE\n"
      " INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE\n";

    //     INITIALIZE
    //     ----------

    if ( iflag == 0 ) goto L10;
    switch ( iflag ) {
    case 1 :
      goto L172;
    case 2 :
      goto L100;
    case 100 :
      goto L_iflag_100_return;
    default :
      throw std::runtime_error("lbfgs: invalid iflag value.");
      break;
    }
  L10:
    iter = 0;
    if ( n <= 0 || m <= 0 ) goto L196;
    if ( gtol <= 1.e-04 ) {
      if ( lp > 0 ) {
        std::printf(
          "\n  GTOL IS LESS THAN OR EQUAL TO 1.e-04\n"
          " IT HAS BEEN RESET TO 9.e-01");
      }
      gtol = 9.e-01;
    }
    nfun = 1;
    point = 0;
    finish = false;
    if ( diagco != 0 ) {
      for ( int i = 1, i_end = n; i <= i_end; ++i ) {
        if ( diag(i) <= zero ) {
          iflag = -2;
          if ( lp > 0 ) std::printf(iflag_minus_2_format, i);
          return;
        }
      }
    }
    else {
      std::fill(diag.begin(), diag.end(), 1.0);
    }

    //     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
    //     ---------------------------------------
    //     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
    //         OTHER TEMPORARY INFORMATION.
    //     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
    //     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
    //         IN THE FORMULA THAT COMPUTES H*G.
    //     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
    //         STEPS.
    //     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
    //         GRADIENT DIFFERENCES.
    //
    //     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
    //     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.

    ispt = n + 2*m;
    iypt = ispt + n*m;
    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      w(ispt+i) = -g(i) * diag(i);
    }
    gnorm = std::sqrt(ddot(n,g,1,g,1));
    if (diagco == 0 || diagco == 1) {
      stp1 = one/gnorm;
    }
    else {
      stp1 = one;
    }

    //     PARAMETERS FOR LINE SEARCH ROUTINE

    ftol = 1.0e-4;
    maxfev = 20;

    if ( iprint(1) >= 0 ) lb1(iprint,iter,nfun,gnorm,n,m,x,f,g,stp,finish);

    //    --------------------
    //     MAIN ITERATION LOOP
    //    --------------------

  L80:
    ++iter;
    info = 0;
    bound = iter - 1;
    if ( iter == 1 ) goto L165;
    if ( iter > m ) bound = m;

    ys = ddot(n,w.get1(iypt+npt+1,n),1,w.get1(ispt+npt+1,n),1);
    if ( diagco == 0 ) {
      double yy = ddot(n,w.get1(iypt+npt+1,n),1,w.get1(iypt+npt+1,n),1);
      for ( int i = 1, i_end = n; i <= i_end; ++i ) {
        diag(i) = ys/yy;
      }
    }
    else {
      iflag = 2;
      return;
    }
  L100:
    if ( diagco != 0 ) {
      for ( int i = 1, i_end = n; i <= i_end; ++i ) {
        if ( diag(i) <= zero ) {
          iflag = -2;
          if ( lp > 0 ) std::printf(iflag_minus_2_format, i);
          return;
        }
      }
    }

    //     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
    //     "Updating quasi-Newton matrices with limited storage",
    //     Mathematics of Computation, Vol.24, No.151, pp. 773-782.f
    //     ---------------------------------------------------------

    cp = point;
    if ( point == 0 ) cp = m;
    w(n+cp) = one/ys;
    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      w(i) = -g(i);
    }
    cp = point;
    for ( int i = 1, i_end = bound; i <= i_end; ++i ) {
      --cp;
      if ( cp == -1 ) cp = m - 1;
      sq = ddot(n,w.get1(ispt + cp*n + 1, n),1,w,1);
      inmc = n + m + cp + 1;
      iycn = iypt + cp*n;
      w(inmc) = w(n + cp + 1)*sq;
      daxpy(n,-w(inmc),w.get1(iycn+1,n),1,w,1);
    }

    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      w(i) *= diag(i);
    }

    for ( int i = 1, i_end = bound; i <= i_end; ++i ) {
      yr = ddot(n,w.get1(iypt + cp*n + 1, n),1,w,1);
      beta = w(n + cp + 1)*yr;
      inmc = n + m + cp + 1;
      beta = w(inmc) - beta;
      iscn = ispt + cp*n;
      daxpy(n,beta,w.get1(iscn + 1, n),1,w,1);
      ++cp;
      if ( cp == m ) cp = 0;
    }

    //     STORE THE NEW SEARCH DIRECTION
    //     ------------------------------

    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      w(ispt + point*n + i) = w(i);
    }

    //     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION
    //     BY USING THE LINE SEARCH ROUTINE MCSRCH
    //     ----------------------------------------------------
  L165:
    nfev = 0;
    stp = one;
    if ( iter == 1 ) stp = stp1;
    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      w(i) = g(i);
    }
    prev_iflag = iflag;
    current_search_direction = w.get1(ispt + point*n + 1, n);
    if (diagco == 2) {
      iflag = 100;
      return;
    }
  L_iflag_100_return:
    iflag = prev_iflag;
    current_search_direction = ref1<double>(0,0);
  L172:
    mcsrch(
      n,x,f,g,w.get1(ispt + point*n + 1, n),
      stp,ftol,xtol,maxfev,info,nfev,diag);
    if ( info == -1 ) {
      iflag = 1;
      return;
    }
    if ( info != 1 ) goto L190;
    nfun += nfev;

    //     COMPUTE THE NEW STEP AND GRADIENT CHANGE
    //     -----------------------------------------

    npt = point*n;
    for ( int i = 1, i_end = n; i <= i_end; ++i ) {
      w(ispt + npt + i) *= stp;
      w(iypt + npt + i) = g(i) - w(i);
    }
    ++point;
    if ( point == m ) point = 0;

    //     TERMINATION TEST
    //     ----------------

    gnorm = std::sqrt(ddot(n,g,1,g,1));
    xnorm = std::sqrt(ddot(n,x,1,x,1));
    xnorm = std::max(1.0,xnorm);
    if ( gnorm/xnorm <= eps ) finish = true;

    if ( iprint(1) >= 0 ) lb1(iprint,iter,nfun,gnorm,n,m,x,f,g,stp,finish);
    if ( finish ) {
      iflag = 0;
      return;
    }
    goto L80;

    //     ------------------------------------------------------------
    //     END OF MAIN ITERATION LOOP. ERROR EXITS.
    //     ------------------------------------------------------------

  L190:
    iflag = -1;
    if ( lp > 0 ) {
      std::printf(
        "\n IFLAG= -1 \n"
        " LINE SEARCH FAILED. SEE DOCUMENTATION OF ROUTINE MCSRCH\n"
        " ERROR RETURN OF LINE SEARCH: INFO= %2d\n"
        " POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT\n"
        " OR INCORRECT TOLERANCES\n", info);
    }
    return;
  L196:
    iflag = -3;
    if ( lp > 0 ) {
      std::printf(
        "\n IFLAG= -3\n"
        " IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)\n");
    }
  }
};

}}} // namespace scitbx::lbfgs::raw

#endif // GUARD
