// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     March 2002: created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_LBFGS_H
#define CCTBX_LBFGS_H

#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <string>
#include <stdio.h>
// The core minimizer can easily be used outside the cctbx and without
// the Boost library if the following includes are removed or
// substituted, depending on the compiler.
#include <boost/config.hpp> // fixes for broken compilers
#include <cctbx/fixes/cmath>
// Includes for drop_convergence_test.
// May be removed if drop_convergence_test is not needed.
#include <cctbx/array_family/shared.h>
#include <cctbx/array_family/tiny_types.h>
#include <cctbx/array_family/tiny_algebra.h>
#include <cctbx/basic/matrixlite.h>
#include <cctbx/math/linear_regression.h>

namespace cctbx {

//! Limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS) %minimizer.
/*! Implementation of the
    Limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS)
    algorithm for large-scale multidimensional minimization
    problems.

    This code was manually derived from Java code which was
    in turn derived from the Fortran program
    <code>lbfgs.f</code>.  The Java translation was
    effected mostly mechanically, with some manual
    clean-up; in particular, array indices start at 0
    instead of 1.  Most of the comments from the Fortran
    code have been pasted in.

    Information on the original LBFGS Fortran source code is
    available at
    http://www.netlib.org/opt/lbfgs_um.shar . The following
    information is taken verbatim from the Netlib documentation
    for the Fortran source.

    <pre>
    file    opt/lbfgs_um.shar
    for     unconstrained optimization problems
    alg     limited memory BFGS method
    by      J. Nocedal
    contact nocedal@eecs.nwu.edu
    ref     D. C. Liu and J. Nocedal, ``On the limited memory BFGS method for
    ,       large scale optimization methods'' Mathematical Programming 45
    ,       (1989), pp. 503-528.
    ,       (Postscript file of this paper is available via anonymous ftp
    ,       to eecs.nwu.edu in the directory pub/%lbfgs/lbfgs_um.)
    </pre>

    @author Jorge Nocedal: original Fortran version, including comments
    (July 1990).<br>
    Robert Dodier: Java translation, August 1997.<br>
    Ralf W. Grosse-Kunstleve: C++ port, March 2002.
 */
namespace lbfgs {

  //! Generic exception class for %lbfgs %error messages.
  /*! All exceptions thrown by the minimizer are derived from this class.
   */
  class error : public std::exception {
    public:
      //! Constructor.
      error(const std::string& msg) throw()
        : msg_("lbfgs error: " + msg)
      {}
      //! Access to error message.
      virtual const char* what() const throw() { return msg_.c_str(); }
    protected:
      virtual ~error() throw() {}
      std::string msg_;
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    public:
      static std::string itoa(unsigned long i) {
        char buf[80];
        sprintf(buf, "%lu", i); // FUTURE: use C++ facility
        return std::string(buf);
      }
#endif
  };

  //! Specific exception class.
  class error_internal_error : public error {
    public:
      //! Constructor.
      error_internal_error(const char* file, unsigned long line) throw()
        : error(
            "Internal Error: " + std::string(file) + "(" + itoa(line) + ")")
      {}
  };

  //! Specific exception class.
  class error_improper_input_parameter : public error {
    public:
      //! Constructor.
      error_improper_input_parameter(const std::string& msg) throw()
        : error("Improper input parameter: " + msg)
      {}
  };

  //! Specific exception class.
  class error_improper_input_data : public error {
    public:
      //! Constructor.
      error_improper_input_data(const std::string& msg) throw()
        : error("Improper input data: " + msg)
      {}
  };

  //! Specific exception class.
  class error_search_direction_not_descent : public error {
    public:
      //! Constructor.
      error_search_direction_not_descent() throw()
        : error("The search direction is not a descent direction.")
      {}
  };

  //! Specific exception class.
  class error_line_search_failed : public error {
    public:
      //! Constructor.
      error_line_search_failed(const std::string& msg) throw()
        : error("Line search failed: " + msg)
      {}
  };

  //! Specific exception class.
  class error_line_search_failed_rounding_errors
  : public error_line_search_failed {
    public:
      //! Constructor.
      error_line_search_failed_rounding_errors(const std::string& msg) throw()
        : error_line_search_failed(msg)
      {}
  };

  namespace detail {

    template <typename NumType>
    inline
    NumType
    pow2(const NumType& x) { return x * x; }

    template <typename NumType>
    inline
    NumType
    abs(const NumType& x) {
      if (x < NumType(0)) return -x;
      return x;
    }

    // This class implements an algorithm for multi-dimensional line search.
    template <typename FloatType>
    class mcsrch
    {
      protected:
        int infoc;
        FloatType dginit;
        bool brackt;
        bool stage1;
        FloatType finit;
        FloatType dgtest;
        FloatType width;
        FloatType width1;
        FloatType stx;
        FloatType fx;
        FloatType dgx;
        FloatType sty;
        FloatType fy;
        FloatType dgy;
        FloatType stmin;
        FloatType stmax;

        static const FloatType& max3(
          const FloatType& x,
          const FloatType& y,
          const FloatType& z)
        {
          return x < y ? (y < z ? z : y ) : (x < z ? z : x );
        }

      public:
        /* Minimize a function along a search direction. This code is
           a Java translation of the function <code>MCSRCH</code> from
           <code>lbfgs.f</code>, which in turn is a slight modification
           of the subroutine <code>CSRCH</code> of More' and Thuente.
           The changes are to allow reverse communication, and do not
           affect the performance of the routine. This function, in turn,
           calls <code>mcstep</code>.<p>

           The Java translation was effected mostly mechanically, with
           some manual clean-up; in particular, array indices start at 0
           instead of 1.  Most of the comments from the Fortran code have
           been pasted in here as well.<p>

           The purpose of <code>mcsrch</code> is to find a step which
           satisfies a sufficient decrease condition and a curvature
           condition.<p>

           At each stage this function updates an interval of uncertainty
           with endpoints <code>stx</code> and <code>sty</code>. The
           interval of uncertainty is initially chosen so that it
           contains a minimizer of the modified function
           <pre>
                f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
           </pre>
           If a step is obtained for which the modified function has a
           nonpositive function value and nonnegative derivative, then
           the interval of uncertainty is chosen so that it contains a
           minimizer of <code>f(x+stp*s)</code>.<p>

           The algorithm is designed to find a step which satisfies
           the sufficient decrease condition
           <pre>
                 f(x+stp*s) &lt;= f(X) + ftol*stp*(gradf(x)'s),
           </pre>
           and the curvature condition
           <pre>
                 abs(gradf(x+stp*s)'s)) &lt;= gtol*abs(gradf(x)'s).
           </pre>
           If <code>ftol</code> is less than <code>gtol</code> and if,
           for example, the function is bounded below, then there is
           always a step which satisfies both conditions. If no step can
           be found which satisfies both conditions, then the algorithm
           usually stops when rounding errors prevent further progress.
           In this case <code>stp</code> only satisfies the sufficient
           decrease condition.<p>

           @author Original Fortran version by Jorge J. More' and
             David J. Thuente as part of the Minpack project, June 1983,
             Argonne National Laboratory. Java translation by Robert
             Dodier, August 1997.

           @param n The number of variables.

           @param x On entry this contains the base point for the line
             search. On exit it contains <code>x + stp*s</code>.

           @param f On entry this contains the value of the objective
             function at <code>x</code>. On exit it contains the value
             of the objective function at <code>x + stp*s</code>.

           @param g On entry this contains the gradient of the objective
             function at <code>x</code>. On exit it contains the gradient
             at <code>x + stp*s</code>.

           @param s The search direction.

           @param stp On entry this contains an initial estimate of a
             satifactory step length. On exit <code>stp</code> contains
             the final estimate.

           @param ftol Tolerance for the sufficient decrease condition.

           @param xtol Termination occurs when the relative width of the
             interval of uncertainty is at most <code>xtol</code>.

           @param maxfev Termination occurs when the number of evaluations
             of the objective function is at least <code>maxfev</code> by
             the end of an iteration.

           @param info This is an output variable, which can have these
             values:
             <ul>
             <li><code>info = -1</code> A return is made to compute
                 the function and gradient.
             <li><code>info = 1</code> The sufficient decrease condition
                 and the directional derivative condition hold.
             </ul>

           @param nfev On exit, this is set to the number of function
             evaluations.

           @param wa Temporary storage array, of length <code>n</code>.
         */
        void run(
          const FloatType& gtol,
          const FloatType& stpmin,
          const FloatType& stpmax,
          std::size_t n,
          FloatType* x,
          FloatType f,
          const FloatType* g,
          FloatType* s,
          std::size_t is0,
          FloatType& stp,
          FloatType ftol,
          FloatType xtol,
          std::size_t maxfev,
          int& info,
          std::size_t& nfev,
          FloatType* wa);

        /* The purpose of this function is to compute a safeguarded step
           for a linesearch and to update an interval of uncertainty for
           a minimizer of the function.<p>

           The parameter <code>stx</code> contains the step with the
           least function value. The parameter <code>stp</code> contains
           the current step. It is assumed that the derivative at
           <code>stx</code> is negative in the direction of the step. If
           <code>brackt</code> is <code>true</code> when
           <code>mcstep</code> returns then a minimizer has been
           bracketed in an interval of uncertainty with endpoints
           <code>stx</code> and <code>sty</code>.<p>

           Variables that must be modified by <code>mcstep</code> are
           implemented as 1-element arrays.

           @param stx Step at the best step obtained so far.
             This variable is modified by <code>mcstep</code>.
           @param fx Function value at the best step obtained so far.
             This variable is modified by <code>mcstep</code>.
           @param dx Derivative at the best step obtained so far.
             The derivative must be negative in the direction of the
             step, that is, <code>dx</code> and <code>stp-stx</code> must
             have opposite signs.  This variable is modified by
             <code>mcstep</code>.

           @param sty Step at the other endpoint of the interval of
             uncertainty. This variable is modified by <code>mcstep</code>.
           @param fy Function value at the other endpoint of the interval
             of uncertainty. This variable is modified by
             <code>mcstep</code>.

           @param dy Derivative at the other endpoint of the interval of
             uncertainty. This variable is modified by <code>mcstep</code>.

           @param stp Step at the current step. If <code>brackt</code> is set
             then on input <code>stp</code> must be between <code>stx</code>
             and <code>sty</code>. On output <code>stp</code> is set to the
             new step.
           @param fp Function value at the current step.
           @param dp Derivative at the current step.

           @param brackt Tells whether a minimizer has been bracketed.
             If the minimizer has not been bracketed, then on input this
             variable must be set <code>false</code>. If the minimizer has
             been bracketed, then on output this variable is
             <code>true</code>.

           @param stpmin Lower bound for the step.
           @param stpmax Upper bound for the step.

           If the return value is 1, 2, 3, or 4, then the step has
           been computed successfully. A return value of 0 indicates
           improper input parameters.

           @author Jorge J. More, David J. Thuente: original Fortran version,
             as part of Minpack project. Argonne Nat'l Laboratory, June 1983.
             Robert Dodier: Java translation, August 1997.
         */
        static int mcstep(
          FloatType& stx,
          FloatType& fx,
          FloatType& dx,
          FloatType& sty,
          FloatType& fy,
          FloatType& dy,
          FloatType& stp,
          FloatType fp,
          FloatType dp,
          bool& brackt,
          FloatType stpmin,
          FloatType stpmax);
    };

    template <typename FloatType>
    void mcsrch<FloatType>::run(
      const FloatType& gtol,
      const FloatType& stpmin,
      const FloatType& stpmax,
      std::size_t n,
      FloatType* x,
      FloatType f,
      const FloatType* g,
      FloatType* s,
      std::size_t is0,
      FloatType& stp,
      FloatType ftol,
      FloatType xtol,
      std::size_t maxfev,
      int& info,
      std::size_t& nfev,
      FloatType* wa)
    {
      if (info != -1) {
        infoc = 1;
        if (   n == 0
            || maxfev == 0
            || gtol < FloatType(0)
            || xtol < FloatType(0)
            || stpmin < FloatType(0)
            || stpmax < stpmin) {
          throw error_internal_error(__FILE__, __LINE__);
        }
        if (stp <= FloatType(0) || ftol < FloatType(0)) {
          throw error_internal_error(__FILE__, __LINE__);
        }
        // Compute the initial gradient in the search direction
        // and check that s is a descent direction.
        dginit = FloatType(0);
        for (std::size_t j = 0; j < n; j++) {
          dginit += g[j] * s[is0+j];
        }
        if (dginit >= FloatType(0)) {
          throw error_search_direction_not_descent();
        }
        brackt = false;
        stage1 = true;
        nfev = 0;
        finit = f;
        dgtest = ftol*dginit;
        width = stpmax - stpmin;
        width1 = FloatType(2) * width;
        std::copy(x, x+n, wa);
        // The variables stx, fx, dgx contain the values of the step,
        // function, and directional derivative at the best step.
        // The variables sty, fy, dgy contain the value of the step,
        // function, and derivative at the other endpoint of
        // the interval of uncertainty.
        // The variables stp, f, dg contain the values of the step,
        // function, and derivative at the current step.
        stx = FloatType(0);
        fx = finit;
        dgx = dginit;
        sty = FloatType(0);
        fy = finit;
        dgy = dginit;
      }
      for (;;) {
        if (info != -1) {
          // Set the minimum and maximum steps to correspond
          // to the present interval of uncertainty.
          if (brackt) {
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
          }
          else {
            stmin = stx;
            stmax = stp + FloatType(4) * (stp - stx);
          }
          // Force the step to be within the bounds stpmax and stpmin.
          stp = std::max(stp, stpmin);
          stp = std::min(stp, stpmax);
          // If an unusual termination is to occur then let
          // stp be the lowest point obtained so far.
          if (   (brackt && (stp <= stmin || stp >= stmax))
              || nfev >= maxfev - 1 || infoc == 0
              || (brackt && stmax - stmin <= xtol * stmax)) {
            stp = stx;
          }
          // Evaluate the function and gradient at stp
          // and compute the directional derivative.
          // We return to main program to obtain F and G.
          for (std::size_t j = 0; j < n; j++) {
            x[j] = wa[j] + stp * s[is0+j];
          }
          info=-1;
          break;
        }
        info = 0;
        nfev++;
        FloatType dg(0);
        for (std::size_t j = 0; j < n; j++) {
          dg += g[j] * s[is0+j];
        }
        FloatType ftest1 = finit + stp*dgtest;
        // Test for convergence.
        if ((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0) {
          throw error_line_search_failed_rounding_errors(
            "Rounding errors prevent further progress."
            " There may not be a step which satisfies the"
            " sufficient decrease and curvature conditions."
            " Tolerances may be too small.");
        }
        if (stp == stpmax && f <= ftest1 && dg <= dgtest) {
          throw error_line_search_failed(
            "The step is at the upper bound stpmax().");
        }
        if (stp == stpmin && (f > ftest1 || dg >= dgtest)) {
          throw error_line_search_failed(
            "The step is at the lower bound stpmin().");
        }
        if (nfev >= maxfev) {
          throw error_line_search_failed(
            "Number of function evaluations has reached maxfev().");
        }
        if (brackt && stmax - stmin <= xtol * stmax) {
          throw error_line_search_failed(
            "Relative width of the interval of uncertainty"
            " is at most xtol().");
        }
        // Check for termination.
        if (f <= ftest1 && abs(dg) <= gtol * (-dginit)) {
          info = 1;
          break;
        }
        // In the first stage we seek a step for which the modified
        // function has a nonpositive value and nonnegative derivative.
        if (   stage1 && f <= ftest1
            && dg >= std::min(ftol, gtol) * dginit) {
          stage1 = false;
        }
        // A modified function is used to predict the step only if
        // we have not obtained a step for which the modified
        // function has a nonpositive function value and nonnegative
        // derivative, and if a lower function value has been
        // obtained but the decrease is not sufficient.
        if (stage1 && f <= fx && f > ftest1) {
          // Define the modified function and derivative values.
          FloatType fm = f - stp*dgtest;
          FloatType fxm = fx - stx*dgtest;
          FloatType fym = fy - sty*dgtest;
          FloatType dgm = dg - dgtest;
          FloatType dgxm = dgx - dgtest;
          FloatType dgym = dgy - dgtest;
          // Call cstep to update the interval of uncertainty
          // and to compute the new step.
          infoc = mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm,
                         brackt, stmin, stmax);
          // Reset the function and gradient values for f.
          fx = fxm + stx*dgtest;
          fy = fym + sty*dgtest;
          dgx = dgxm + dgtest;
          dgy = dgym + dgtest;
        }
        else {
          // Call mcstep to update the interval of uncertainty
          // and to compute the new step.
          infoc = mcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg,
                         brackt, stmin, stmax);
        }
        // Force a sufficient decrease in the size of the
        // interval of uncertainty.
        if (brackt) {
          if (abs(sty - stx) >= FloatType(0.66) * width1) {
            stp = stx + FloatType(0.5) * (sty - stx);
          }
          width1 = width;
          width = abs(sty - stx);
        }
      }
    }

    template <typename FloatType>
    int mcsrch<FloatType>::mcstep(
      FloatType& stx,
      FloatType& fx,
      FloatType& dx,
      FloatType& sty,
      FloatType& fy,
      FloatType& dy,
      FloatType& stp,
      FloatType fp,
      FloatType dp,
      bool& brackt,
      FloatType stpmin,
      FloatType stpmax)
    {
      bool bound;
      FloatType gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
      int info = 0;
      if (   (   brackt && (stp <= std::min(stx, sty)
              || stp >= std::max(stx, sty)))
          || dx * (stp - stx) >= FloatType(0) || stpmax < stpmin) {
        return 0;
      }
      // Determine if the derivatives have opposite sign.
      sgnd = dp * (dx / abs(dx));
      if (fp > fx) {
        // First case. A higher function value.
        // The minimum is bracketed. If the cubic step is closer
        // to stx than the quadratic step, the cubic step is taken,
        // else the average of the cubic and quadratic steps is taken.
        info = 1;
        bound = true;
        theta = FloatType(3) * (fx - fp) / (stp - stx) + dx + dp;
        s = max3(abs(theta), abs(dx), abs(dp));
        gamma = s * std::sqrt(pow2(theta / s) - (dx / s) * (dp / s));
        if (stp < stx) gamma = - gamma;
        p = (gamma - dx) + theta;
        q = ((gamma - dx) + gamma) + dp;
        r = p/q;
        stpc = stx + r * (stp - stx);
        stpq = stx
          + ((dx / ((fx - fp) / (stp - stx) + dx)) / FloatType(2))
            * (stp - stx);
        if (abs(stpc - stx) < abs(stpq - stx)) {
          stpf = stpc;
        }
        else {
          stpf = stpc + (stpq - stpc) / FloatType(2);
        }
        brackt = true;
      }
      else if (sgnd < FloatType(0)) {
        // Second case. A lower function value and derivatives of
        // opposite sign. The minimum is bracketed. If the cubic
        // step is closer to stx than the quadratic (secant) step,
        // the cubic step is taken, else the quadratic step is taken.
        info = 2;
        bound = false;
        theta = FloatType(3) * (fx - fp) / (stp - stx) + dx + dp;
        s = max3(abs(theta), abs(dx), abs(dp));
        gamma = s * std::sqrt(pow2(theta / s) - (dx / s) * (dp / s));
        if (stp > stx) gamma = - gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dx;
        r = p/q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (abs(stpc - stp) > abs(stpq - stp)) {
          stpf = stpc;
        }
        else {
          stpf = stpq;
        }
        brackt = true;
      }
      else if (abs(dp) < abs(dx)) {
        // Third case. A lower function value, derivatives of the
        // same sign, and the magnitude of the derivative decreases.
        // The cubic step is only used if the cubic tends to infinity
        // in the direction of the step or if the minimum of the cubic
        // is beyond stp. Otherwise the cubic step is defined to be
        // either stpmin or stpmax. The quadratic (secant) step is also
        // computed and if the minimum is bracketed then the the step
        // closest to stx is taken, else the step farthest away is taken.
        info = 3;
        bound = true;
        theta = FloatType(3) * (fx - fp) / (stp - stx) + dx + dp;
        s = max3(abs(theta), abs(dx), abs(dp));
        gamma = s * std::sqrt(
          std::max(FloatType(0), pow2(theta / s) - (dx / s) * (dp / s)));
        if (stp > stx) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = (gamma + (dx - dp)) + gamma;
        r = p/q;
        if (r < FloatType(0) && gamma != FloatType(0)) {
          stpc = stp + r * (stx - stp);
        }
        else if (stp > stx) {
          stpc = stpmax;
        }
        else {
          stpc = stpmin;
        }
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (brackt) {
          if (abs(stp - stpc) < abs(stp - stpq)) {
            stpf = stpc;
          }
          else {
            stpf = stpq;
          }
        }
        else {
          if (abs(stp - stpc) > abs(stp - stpq)) {
            stpf = stpc;
          }
          else {
            stpf = stpq;
          }
        }
      }
      else {
        // Fourth case. A lower function value, derivatives of the
        // same sign, and the magnitude of the derivative does
        // not decrease. If the minimum is not bracketed, the step
        // is either stpmin or stpmax, else the cubic step is taken.
        info = 4;
        bound = false;
        if (brackt) {
          theta = FloatType(3) * (fp - fy) / (sty - stp) + dy + dp;
          s = max3(abs(theta), abs(dy), abs(dp));
          gamma = s * std::sqrt(pow2(theta / s) - (dy / s) * (dp / s));
          if (stp > sty) gamma = -gamma;
          p = (gamma - dp) + theta;
          q = ((gamma - dp) + gamma) + dy;
          r = p/q;
          stpc = stp + r * (sty - stp);
          stpf = stpc;
        }
        else if (stp > stx) {
          stpf = stpmax;
        }
        else {
          stpf = stpmin;
        }
      }
      // Update the interval of uncertainty. This update does not
      // depend on the new step or the case analysis above.
      if (fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
      }
      else {
        if (sgnd < FloatType(0)) {
          sty = stx;
          fy = fx;
          dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
      }
      // Compute the new step and safeguard it.
      stpf = std::min(stpmax, stpf);
      stpf = std::max(stpmin, stpf);
      stp = stpf;
      if (brackt && bound) {
        if (sty > stx) {
          stp = std::min(stx + FloatType(0.66) * (sty - stx), stp);
        }
        else {
          stp = std::max(stx + FloatType(0.66) * (sty - stx), stp);
        }
      }
      return info;
    }

    /* Compute the sum of a vector times a scalar plus another vector.
       Adapted from the subroutine <code>daxpy</code> in
       <code>lbfgs.f</code>.
     */
    template <typename FloatType>
    void daxpy(
      std::size_t n,
      FloatType da,
      const FloatType* dx,
      std::size_t ix0,
      std::size_t incx,
      FloatType* dy,
      std::size_t iy0,
      std::size_t incy)
    {
      std::size_t i, ix, iy, m;
      if (n == 0) return;
      if (da == FloatType(0)) return;
      if  (!(incx == 1 && incy == 1)) {
        ix = 0;
        iy = 0;
        for (i = 0; i < n; i++) {
          dy[iy0+iy] += da * dx[ix0+ix];
          ix += incx;
          iy += incy;
        }
        return;
      }
      m = n % 4;
      for (i = 0; i < m; i++) {
        dy[iy0+i] += da * dx[ix0+i];
      }
      for (; i < n;) {
        dy[iy0+i] += da * dx[ix0+i]; i++;
        dy[iy0+i] += da * dx[ix0+i]; i++;
        dy[iy0+i] += da * dx[ix0+i]; i++;
        dy[iy0+i] += da * dx[ix0+i]; i++;
      }
    }

    /* Compute the dot product of two vectors.
       Adapted from the subroutine <code>ddot</code>
       in <code>lbfgs.f</code>.
     */
    template <typename FloatType>
    FloatType ddot(
      std::size_t n,
      const FloatType* dx,
      std::size_t ix0,
      std::size_t incx,
      const FloatType* dy,
      std::size_t iy0,
      std::size_t incy)
    {
      std::size_t i, ix, iy, m, mp1;
      FloatType dtemp(0);
      if (n == 0) return FloatType(0);
      if (!(incx == 1 && incy == 1)) {
        ix = 0;
        iy = 0;
        for (i = 0; i < n; i++) {
          dtemp += dx[ix0+ix] * dy[iy0+iy];
          ix += incx;
          iy += incy;
        }
        return dtemp;
      }
      m = n % 5;
      for (i = 0; i < m; i++) {
        dtemp += dx[ix0+i] * dy[iy0+i];
      }
      for (; i < n;) {
        dtemp += dx[ix0+i] * dy[iy0+i]; i++;
        dtemp += dx[ix0+i] * dy[iy0+i]; i++;
        dtemp += dx[ix0+i] * dy[iy0+i]; i++;
        dtemp += dx[ix0+i] * dy[iy0+i]; i++;
        dtemp += dx[ix0+i] * dy[iy0+i]; i++;
      }
      return dtemp;
    }

  } // namespace detail

  //! Interface to the LBFGS %minimizer.
  /*! This class solves the unconstrained minimization problem
      <pre>
          min f(x),  x = (x1,x2,...,x_n),
      </pre>
      using the limited-memory BFGS method. The routine is
      especially effective on problems involving a large number of
      variables. In a typical iteration of this method an
      approximation Hk to the inverse of the Hessian
      is obtained by applying <code>m</code> BFGS updates to a
      diagonal matrix Hk0, using information from the
      previous <code>m</code> steps.  The user specifies the number
      <code>m</code>, which determines the amount of storage
      required by the routine. The user may also provide the
      diagonal matrices Hk0 (parameter <code>diag</code> in the run()
      function) if not satisfied with the default choice. The
      algorithm is described in "On the limited memory BFGS method for
      large scale optimization", by D. Liu and J. Nocedal, Mathematical
      Programming B 45 (1989) 503-528.

      The user is required to calculate the function value
      <code>f</code> and its gradient <code>g</code>. In order to
      allow the user complete control over these computations,
      reverse communication is used. The routine must be called
      repeatedly under the control of the member functions
      <code>requests_f_and_g()</code>,
      <code>requests_diag()</code>.
      If neither requests_f_and_g() nor requests_diag() is
      <code>true</code> the user should check for convergence
      (using class traditional_convergence_test or any
      other custom test). If the convergence test is negative,
      the minimizer may be called again for the next iteration.

      The steplength (stp()) is determined at each iteration
      by means of the line search routine <code>mcsrch</code>, which is
      a slight modification of the routine <code>CSRCH</code> written
      by More' and Thuente.

      The only variables that are machine-dependent are
      <code>xtol</code>,
      <code>stpmin</code> and
      <code>stpmax</code>.

      Fatal errors cause <code>error</code> exceptions to be thrown.
      The generic class <code>error</code> is sub-classed (e.g.
      class <code>error_line_search_failed</code>) to facilitate
      granular %error handling.

      A note on performance: Using Compaq Fortran V5.4 and
      Compaq C++ V6.5, the C++ implementation is about 15% slower
      than the Fortran implementation.
   */
  template <typename FloatType>
  class minimizer
  {
    public:
      //! Default constructor. Some members are not initialized!
      minimizer()
      : n_(0), m_(0), maxfev_(0),
        gtol_(0), xtol_(0),
        stpmin_(0), stpmax_(0),
        ispt(0), iypt(0)
      {}

      //! Constructor.
      /*! @param n The number of variables in the minimization problem.
             Restriction: <code>n &gt; 0</code>.

          @param m The number of corrections used in the BFGS update.
             Values of <code>m</code> less than 3 are not recommended;
             large values of <code>m</code> will result in excessive
             computing time. <code>3 &lt;= m &lt;= 7</code> is
             recommended.
             Restriction: <code>m &gt; 0</code>.

          @param maxfev Termination occurs when the number of evaluations
             of the objective function is at least <code>maxfev</code> by
             the end of an iteration.

          @param gtol Controls the accuracy of the line search.
            If the function and gradient evaluations are inexpensive with
            respect to the cost of the iteration (which is sometimes the
            case when solving very large problems) it may be advantageous
            to set <code>gtol</code> to a small value. A typical small
            value is 0.1.
            Restriction: <code>gtol</code> should be greater than 1e-4.

          @param xtol An estimate of the machine precision (e.g. 10e-16
            on a SUN station 3/60). The line search routine will
            terminate if the relative width of the interval of
            uncertainty is less than <code>xtol</code>.

          @param stpmin Specifies the lower bound for the step
            in the line search.
            The default value is 1e-20. This value need not be modified
            unless the exponent is too large for the machine being used,
            or unless the problem is extremely badly scaled (in which
            case the exponent should be increased).

          @param stpmax specifies the upper bound for the step
            in the line search.
            The default value is 1e20. This value need not be modified
            unless the exponent is too large for the machine being used,
            or unless the problem is extremely badly scaled (in which
            case the exponent should be increased).
       */
      explicit
      minimizer(
        std::size_t n,
        std::size_t m = 5,
        std::size_t maxfev = 20,
        FloatType gtol = FloatType(0.9),
        FloatType xtol = FloatType(1.e-16),
        FloatType stpmin = FloatType(1.e-20),
        FloatType stpmax = FloatType(1.e20))
        : n_(n), m_(m), maxfev_(maxfev),
          gtol_(gtol), xtol_(xtol),
          stpmin_(stpmin), stpmax_(stpmax),
          iflag_(0), requests_f_and_g_(false), requests_diag_(false),
          iter_(0), nfun_(0), stp_(0),
          stp1(0), ftol(0.0001), ys(0), point(0), npt(0),
          ispt(n+2*m), iypt((n+2*m)+n*m),
          info(0), bound(0), nfev(0)
      {
        if (n_ == 0) {
          throw error_improper_input_parameter("n = 0.");
        }
        if (m_ == 0) {
          throw error_improper_input_parameter("m = 0.");
        }
        if (maxfev_ == 0) {
         throw error_improper_input_parameter("maxfev = 0.");
        }
        if (gtol_ <= FloatType(1.e-4)) {
          throw error_improper_input_parameter("gtol <= 1.e-4.");
        }
        if (xtol_ < FloatType(0)) {
          throw error_improper_input_parameter("xtol < 0.");
        }
        if (stpmin_ < FloatType(0)) {
          throw error_improper_input_parameter("stpmin < 0.");
        }
        if (stpmax_ < stpmin) {
          throw error_improper_input_parameter("stpmax < stpmin");
        }
        w_.resize(n_*(2*m_+1)+2*m_);
      }

      //! Number of free parameters (as passed to the constructor).
      std::size_t n() const { return n_; }

      //! Number of corrections kept (as passed to the constructor).
      std::size_t m() const { return m_; }

      /*! \brief Maximum number of evaluations of the objective function
          (as passed to the constructor).
       */
      std::size_t maxfev() const { return maxfev_; }

      /*! \brief Control of the accuracy of the line search.
          (as passed to the constructor).
       */
      FloatType gtol() const { return gtol_; }

      //! Estimate of the machine precision (as passed to the constructor).
      FloatType xtol() const { return xtol_; }

      /*! \brief Lower bound for the step in the line search.
          (as passed to the constructor).
       */
      FloatType stpmin() const { return stpmin_; }

      /*! \brief Upper bound for the step in the line search.
          (as passed to the constructor).
       */
      FloatType stpmax() const { return stpmax_; }

      //! Status indicator for reverse communication.
      /*! <code>true</code> if the run() function returns to request
          evaluation of the objective function (<code>f</code>) and
          gradients (<code>g</code>) for the current point
          (<code>x</code>). To continue the minimization the
          run() function is called again with the updated values for
          <code>f</code> and <code>g</code>.
          <p>
          See also: requests_diag()
       */
      bool requests_f_and_g() const { return requests_f_and_g_; }

      //! Status indicator for reverse communication.
      /*! <code>true</code> if the run() function returns to request
          evaluation of the diagonal matrix (<code>diag</code>)
          for the current point (<code>x</code>).
          To continue the minimization the run() function is called
          again with the updated values for <code>diag</code>.
          <p>
          See also: requests_f_and_g()
       */
      bool requests_diag() const { return requests_diag_; }

      //! Number of iterations so far.
      /*! Note that one iteration may involve multiple evaluations
          of the objective function.
          <p>
          See also: nfun()
       */
      std::size_t iter() const { return iter_; }

      //! Total number of evaluations of the objective function so far.
      /*! The total number of function evaluations increases by the
          number of evaluations required for the line search. The total
          is only increased after a successful line search.
          <p>
          See also: iter()
       */
      std::size_t nfun() const { return nfun_; }

      //! Norm of gradient given gradient array of length n().
      FloatType euclidean_norm(const FloatType* a) const {
        return std::sqrt(detail::ddot(n_, a, 0, 1, a, 0, 1));
      }

      //! Current stepsize.
      FloatType stp() const { return stp_; }

      //! Execution of one step of the minimization.
      /*! @param x On initial entry this must be set by the user to
             the values of the initial estimate of the solution vector.

          @param f Before initial entry or on re-entry under the
             control of requests_f_and_g(), <code>f</code> must be set
             by the user to contain the value of the objective function
             at the current point <code>x</code>.

          @param g Before initial entry or on re-entry under the
             control of requests_f_and_g(), <code>g</code> must be set
             by the user to contain the components of the gradient at
             the current point <code>x</code>.

          The return value is <code>true</code> if either
          requests_f_and_g() or requests_diag() is <code>true</code>.
          Otherwise the user should check for convergence
          (e.g. using class traditional_convergence_test) and
          call the run() function again to continue the minimization.
          If the return value is <code>false</code> the user
          should <b>not</b> update <code>f</code>, <code>g</code> or
          <code>diag</code> (other overload) before calling
          the run() function again.

          Note that <code>x</code> is always modified by the run()
          function. Depending on the situation it can therefore be
          necessary to evaluate the objective function one more time
          after the minimization is terminated.
       */
      bool run(
        FloatType* x,
        FloatType f,
        const FloatType* g)
      {
        if (diag_.size() == 0) diag_.resize(n_);
        return generic_run(x, f, g, false, &(*(diag_.begin())));
      }

      //! Execution of one step of the minimization.
      /*! @param x See other overload.

          @param f See other overload.

          @param g See other overload.

          @param diag On initial entry or on re-entry under the
             control of requests_diag(), <code>diag</code> must be set by
             the user to contain the values of the diagonal matrix Hk0.
             The routine will return at each iteration of the algorithm
             with requests_diag() set to <code>true</code>.
             <p>
             Restriction: all elements of <code>diag</code> must be
             positive.
       */
      bool run(
        FloatType* x,
        FloatType f,
        const FloatType* g,
        FloatType* diag)
      {
        return generic_run(x, f, g, true, diag);
      }

    protected:
      static void throw_diagonal_element_not_positive(std::size_t i) {
        throw error_improper_input_data(
          "The " + error::itoa(i) + ". diagonal element of the"
          " inverse Hessian approximation is not positive.");
      }

      bool generic_run(
        FloatType* x,
        FloatType f,
        const FloatType* g,
        bool diagco,
        FloatType* diag);

      detail::mcsrch<FloatType> mcsrch_instance;
      const std::size_t n_;
      const std::size_t m_;
      const std::size_t maxfev_;
      const FloatType gtol_;
      const FloatType xtol_;
      const FloatType stpmin_;
      const FloatType stpmax_;
      int iflag_;
      bool requests_f_and_g_;
      bool requests_diag_;
      std::size_t iter_;
      std::size_t nfun_;
      FloatType stp_;
      FloatType stp1;
      FloatType ftol;
      FloatType ys;
      std::size_t point;
      std::size_t npt;
      const std::size_t ispt;
      const std::size_t iypt;
      int info;
      std::size_t bound;
      std::size_t nfev;
      std::vector<FloatType> w_;
      std::vector<FloatType> diag_;
  };

  template <typename FloatType>
  bool minimizer<FloatType>::generic_run(
    FloatType* x,
    FloatType f,
    const FloatType* g,
    bool diagco,
    FloatType* diag)
  {
    bool execute_entire_while_loop = false;
    if (!(requests_f_and_g_ || requests_diag_)) {
      execute_entire_while_loop = true;
    }
    requests_f_and_g_ = false;
    requests_diag_ = false;
    FloatType* w = &(*(w_.begin()));
    if (iflag_ == 0) { // Initialize.
      nfun_ = 1;
      if (diagco) {
        for (std::size_t i = 0; i < n_; i++) {
          if (diag[i] <= FloatType(0)) {
            throw_diagonal_element_not_positive(i);
          }
        }
      }
      else {
        std::fill_n(diag, n_, FloatType(1));
      }
      for (std::size_t i = 0; i < n_; i++) {
        w[ispt + i] = -g[i] * diag[i];
      }
      FloatType gnorm = std::sqrt(detail::ddot(n_, g, 0, 1, g, 0, 1));
      if (gnorm == FloatType(0)) return false;
      stp1 = FloatType(1) / gnorm;
      execute_entire_while_loop = true;
    }
    if (execute_entire_while_loop) {
      bound = iter_;
      iter_++;
      info = 0;
      if (iter_ != 1) {
        if (iter_ > m_) bound = m_;
        ys = detail::ddot(n_, w, iypt + npt, 1, w, ispt + npt, 1);
        if (!diagco) {
          FloatType yy = detail::ddot(n_, w, iypt+npt, 1, w, iypt+npt, 1);
          std::fill_n(diag, n_, ys / yy);
        }
        else {
          iflag_ = 2;
          requests_diag_ = true;
          return true;
        }
      }
    }
    if (execute_entire_while_loop || iflag_ == 2) {
      if (iter_ != 1) {
        if (diagco) {
          for (std::size_t i = 0; i < n_; i++) {
            if (diag[i] <= FloatType(0)) {
              throw_diagonal_element_not_positive(i);
            }
          }
        }
        std::size_t cp = point;
        if (point == 0) cp = m_;
        w[n_ + cp -1] = 1 / ys;
        std::size_t i;
        for (i = 0; i < n_; i++) {
          w[i] = -g[i];
        }
        cp = point;
        for (i = 0; i < bound; i++) {
          if (cp == 0) cp = m_;
          cp--;
          FloatType sq = detail::ddot(n_, w, ispt + cp * n_, 1, w, 0, 1);
          std::size_t inmc=n_+m_+cp;
          std::size_t iycn=iypt+cp*n_;
          w[inmc] = w[n_ + cp] * sq;
          detail::daxpy(n_, -w[inmc], w, iycn, 1, w, 0, 1);
        }
        for (i = 0; i < n_; i++) {
          w[i] *= diag[i];
        }
        for (i = 0; i < bound; i++) {
          FloatType yr = detail::ddot(n_, w, iypt + cp * n_, 1, w, 0, 1);
          FloatType beta = w[n_ + cp] * yr;
          std::size_t inmc=n_+m_+cp;
          beta = w[inmc] - beta;
          std::size_t iscn=ispt+cp*n_;
          detail::daxpy(n_, beta, w, iscn, 1, w, 0, 1);
          cp++;
          if (cp == m_) cp = 0;
        }
        std::copy(w, w+n_, w+(ispt + point * n_));
      }
      stp_ = FloatType(1);
      if (iter_ == 1) stp_ = stp1;
      std::copy(g, g+n_, w);
    }
    mcsrch_instance.run(
      gtol_, stpmin_, stpmax_, n_, x, f, g, w, ispt + point * n_,
      stp_, ftol, xtol_, maxfev_, info, nfev, diag);
    if (info == -1) {
      iflag_ = 1;
      requests_f_and_g_ = true;
      return true;
    }
    if (info != 1) {
      throw error_internal_error(__FILE__, __LINE__);
    }
    nfun_ += nfev;
    npt = point*n_;
    for (std::size_t i = 0; i < n_; i++) {
      w[ispt + npt + i] = stp_ * w[ispt + npt + i];
      w[iypt + npt + i] = g[i] - w[i];
    }
    point++;
    if (point == m_) point = 0;
    return false;
  }

  //! Traditional LBFGS convergence test.
  /*! This convergence test is equivalent to the test embedded
      in the <code>lbfgs.f</code> Fortran code. The test assumes that
      there is a meaningful relation between the Euclidean norm of the
      parameter vector <code>x</code> and the norm of the gradient
      vector <code>g</code>. Therefore this test should not be used if
      this assumption is not correct for a given problem.
   */
  template <typename FloatType>
  class traditional_convergence_test
  {
    public:
      //! Default constructor.
      traditional_convergence_test()
      : n_(0), eps_(0)
      {}

      //! Constructor.
      /*! @param n The number of variables in the minimization problem.
             Restriction: <code>n &gt; 0</code>.

          @param eps Determines the accuracy with which the solution
            is to be found.
       */
      explicit
      traditional_convergence_test(
        std::size_t n,
        FloatType eps = FloatType(1.e-5))
      : n_(n), eps_(eps)
      {
        if (n_ == 0) {
          throw error_improper_input_parameter("n = 0.");
        }
        if (eps_ < FloatType(0)) {
          throw error_improper_input_parameter("eps < 0.");
        }
      }

      //! Number of free parameters (as passed to the constructor).
      std::size_t n() const { return n_; }

      /*! \brief Accuracy with which the solution is to be found
          (as passed to the constructor).
       */
      FloatType eps() const { return eps_; }

      //! Execution of the convergence test for the given parameters.
      /*! Returns <code>true</code> if
          <pre>
            ||g|| &lt; eps * max(1,||x||),
          </pre>
          where <code>||.||</code> denotes the Euclidean norm.

          @param x Current solution vector.

          @param g Components of the gradient at the current
            point <code>x</code>.
       */
      bool
      operator()(const FloatType* x, const FloatType* g) const
      {
        FloatType xnorm = std::sqrt(detail::ddot(n_, x, 0, 1, x, 0, 1));
        FloatType gnorm = std::sqrt(detail::ddot(n_, g, 0, 1, g, 0, 1));
        if (gnorm <= eps_ * std::max(FloatType(1), xnorm)) return true;
        return false;
      }
    protected:
      const std::size_t n_;
      const FloatType eps_;
  };

  //! Convergence test based on monitoring the objective function.
  /*! This test monitors the behavior of objective function <code>f</code>
      in the course of the minimization. A combination of criteria
      is used to determine when the objective function has
      reached a plateau without making assumptions about
      the absolute values of the objective function:

      <ul>
      <li>The maximum value of the drop in the objective function
          is determined (max_drop()).
          <p>
      <li>A straight line is fitted to the last p() values
          of the objective function. Let <code>slope1</code>
          be the slope of this fitted line. A necessary condition
          for the detection of convergence is:
          <pre>
          -slope &lt; max_drop * max_drop_eps * number_of_iterations^2
          </pre>
          Note that this test becomes increasingly
          tolerant with the number of iterations. The rationale
          is that after a large number of iterations the second test
          desribed below should be dominant.
          <p>
      <li>If the test for the slope passes, all points of the
          objective function are fitted to
          an exponential curve: f = r * exp(s * number_of_iteration).
          <p>
      <li>The exponential fit is used to compute idealized
          values for the objective function at the
          last p() iterations. Let <code>slope2</code> be
          the slope of this fitted line.
          <p>
      <li>A straight line is fitted to the idealized values.
          Convergence is detected if the following is true:
          <pre>
          abs(slope2 - slope1) < slope_eps
          </pre>
      </ul>
   */
  template <typename FloatType>
  class drop_convergence_test
  {
    public:
      //! Constructor.
      /*! See the class details for an outline of the convergence
          detection algorithm and the meaning of the parameters.
       */
      explicit
      drop_convergence_test(
        std::size_t p = 5,
        FloatType max_drop_eps = FloatType(1.e-5),
        FloatType slope_eps = FloatType(1.e-4))
      : p_(p), max_drop_eps_(max_drop_eps), slope_eps_(slope_eps),
        max_drop_(0)
      {
        if (p_ < 2) {
          throw error_improper_input_parameter("p < 2.");
        }
        if (max_drop_eps_ < FloatType(0)) {
          throw error_improper_input_parameter("max_drop_eps < 0.");
        }
        if (slope_eps_ < FloatType(0)) {
          throw error_improper_input_parameter("slope_eps < 0.");
        }
      }

      /*! \brief Number of most recent objective function values used
           in fit of straight lines (as passed to the constructor).
       */
      std::size_t p() const { return p_; }

      /*! \brief Base tolerance for test of drop of objective function
           (as passed to the constructor).
       */
      FloatType max_drop_eps() const { return max_drop_eps_; }

      /*! \brief Tolerance for comparison of slopes
           (as passed to the constructor).
       */
      FloatType slope_eps() const { return slope_eps_; }

      //! Execution of the convergence test.
      /*! Note that at least p() executions are required before
          convergence will be evaluated.
       */
      bool
      operator()(FloatType f);

      //! Values of the objective function monitored so far.
      af::shared<FloatType> objective_function_values() const {
        return y_;
      }

      //! Maximum drop monitored so far.
      FloatType max_drop() const { return max_drop_; }

    protected:
      const std::size_t p_;
      const FloatType max_drop_eps_;
      const FloatType slope_eps_;
      af::shared<FloatType> x_;
      af::shared<FloatType> y_;
      af::shared<FloatType> ln_y_;
      FloatType max_drop_;
  };

  template <typename FloatType>
  bool
  drop_convergence_test<FloatType>::operator()(FloatType f)
  {
    if (x_.size()) {
      max_drop_ = std::max(max_drop_, y_.back() - f);
    }
    x_.push_back(x_.size() + 1);
    y_.push_back(f);
    ln_y_.push_back(std::log(f));
    if (x_.size() < p_) return false;
    // fit last p_ points to straight line: y = m*x + b
    math::linear_regression<FloatType> linreg_y(
      af::const_ref<FloatType>(x_.end() - p_, p_),
      af::const_ref<FloatType>(y_.end() - p_, p_));
    cctbx_assert(linreg_y.is_well_defined());
    // check absolute value of slope
    FloatType sliding_tolerance =
      max_drop_ * max_drop_eps_ * detail::pow2(FloatType(x_.size()));
    if (-linreg_y.m() > sliding_tolerance) {
      return false;
    }
    // fit all points to exponential: y = r * exp(s * x)
    math::linear_regression<FloatType> linreg_ln_y(
      x_.const_ref(), ln_y_.const_ref());
    cctbx_assert(linreg_ln_y.is_well_defined());
    FloatType r = std::exp(linreg_ln_y.b());
    FloatType s = linreg_ln_y.m();
    // compute exponentially fitted values for last p_ points
    af::shared<FloatType> ideal_y;
    for(std::size_t i = x_.size() - p_; i < x_.size(); i++) {
      ideal_y.push_back(r * std::exp(s * x_[i]));
    }
    // fit computed values to straight line
    math::linear_regression<FloatType> linreg_ideal_y(
      af::const_ref<FloatType>(x_.end() - p_, p_),
      ideal_y.const_ref());
    cctbx_assert(linreg_ideal_y.is_well_defined());
    // compare slopes
    if (detail::abs(linreg_ideal_y.m() - linreg_y.m()) < slope_eps_) {
      return true;
    }
    return false;
  }

}} // namespace cctbx::lbfgs

#endif // CCTBX_LBFGS_H
