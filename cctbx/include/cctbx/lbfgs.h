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

/*! \file
    This file contains code for the
    Limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS)
    algorithm for large-scale multidimensional minimization
    problems.

    This code was manually derived from Java code which was
    in turn derived from the Fortran program
    <code>lbfgs.f</code>.  The Java translation was
    effected mostly mechanically, with some manual
    clean-up; in particular, array indices start at 0
    instead of 1.  Most of the comments from the Fortran
    code have been pasted in here as well.

    <p>
    Information on the original LBFGS Fortran source code is
    available at
    http://www.netlib.org/opt/lbfgs_um.shar . The following
    information is taken verbatim from the Netlib documentation
    for the Fortran source.

    <p>
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
    ,       to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_um.)
    </pre>

    @author Jorge Nocedal: original Fortran version, including comments
    (July 1990). Robert Dodier: Java translation, August 1997.
    Ralf W. Grosse-Kunstleve, C++ port, March 2002.
 */

#include <vector>
#include <string>
#include <stdio.h>
#include <boost/config.hpp>
#include <cctbx/error.h>
#include <cctbx/fixes/cstdlib>
#include <cctbx/fixes/cmath>
#include <cctbx/array_family/shared.h>

namespace cctbx {

  std::string itoa(int i) {
    char buf[80];
    sprintf(buf, "%d", i); // FUTURE: use C++ facility
    return std::string(buf);
  }

  std::string ftoa(double f) {
    char buf[80];
    sprintf(buf, "%10.4g", f); // FUTURE: use C++ facility
    return std::string(buf);
  }

  std::string operator+(const std::string& s, int i) {
    return s + itoa(i);
  }

  std::string operator+(const std::string& s, double f) {
    return s + ftoa(f);
  }

  class lbfgs_parameters {
    public:
      explicit
      lbfgs_parameters(
        double gtol_ = 0.9,
        double stpmin_ = 1.e-20,
        double stpmax_ = 1.e20)
        : gtol(gtol_),
          stpmin(stpmin_),
          stpmax(stpmax_)
      {}
      double gtol;
      double stpmin;
      double stpmax;
  };

  // This class implements an algorithm for multi-dimensional line search.
  class mcsrch
  {
    protected:
      int infoc[1];
      int j;
      double dg;
      double dgm;
      double dginit;
      double dgtest;
      double dgx[1];
      double dgxm[1];
      double dgy[1];
      double dgym[1];
      double finit;
      double ftest1;
      double fm;
      double fx[1];
      double fxm[1];
      double fy[1];
      double fym[1];
      double p5;
      double p66;
      double stx[1];
      double sty[1];
      double stmin;
      double stmax;
      double width;
      double width1;
      double xtrapf;
      bool brackt[1];
      bool stage1;

      void init()
      {
        infoc[0] = 0;
        j = 0;
        dg = 0;
        dgm = 0;
        dginit = 0;
        dgtest = 0;
        dgx[0] = 0;
        dgxm[0] = 0;
        dgy[0] = 0;
        dgym[0] = 0;
        finit = 0;
        ftest1 = 0;
        fm = 0;
        fx[0] = 0;
        fxm[0] = 0;
        fy[0] = 0;
        fym[0] = 0;
        p5 = 0;
        p66 = 0;
        stx[0] = 0;
        sty[0] = 0;
        stmin = 0;
        stmax = 0;
        width = 0;
        width1 = 0;
        xtrapf = 0;
        brackt[0] = 0;
        stage1 = false;
      }

    public:
      static double sqr(double x) { return x * x; }

      static double max3(double x, double y, double z) {
        return x < y ? (y < z ? z : y ) : (x < z ? z : x );
      }

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

         @param x On entry this contains the base point for the line search.
           On exit it contains <code>x + stp*s</code>.

         @param f On entry this contains the value of the objective function
           at <code>x</code>. On exit it contains the value of the objective
           function at <code>x + stp*s</code>.

         @param g On entry this contains the gradient of the objective function
           at <code>x</code>. On exit it contains the gradient at
           <code>x + stp*s</code>.

         @param s The search direction.

         @param stp On entry this contains an initial estimate of a satifactory
           step length. On exit <code>stp</code> contains the final estimate.

         @param ftol Tolerance for the sufficient decrease condition.

         @param xtol Termination occurs when the relative width of the interval
           of uncertainty is at most <code>xtol</code>.

         @param maxfev Termination occurs when the number of evaluations of
           the objective function is at least <code>maxfev</code> by the end
           of an iteration.

         @param info This is an output variable, which can have these values:
           <ul>
           <li><code>info = 0</code> Improper input parameters.
           <li><code>info = -1</code> A return is made to compute
               the function and gradient.
           <li><code>info = 1</code> The sufficient decrease condition
               and the directional derivative condition hold.
           <li><code>info = 2</code> Relative width of the interval of
               uncertainty is at most <code>xtol</code>.
           <li><code>info = 3</code> Number of function evaluations
               has reached <code>maxfev</code>.
           <li><code>info = 4</code> The step is at the lower bound
               <code>stpmin</code>.
           <li><code>info = 5</code> The step is at the upper bound
               <code>stpmax</code>.
           <li><code>info = 6</code> Rounding errors prevent further
               progress. There may not be a step which satisfies the
               sufficient decrease and curvature conditions.
               Tolerances may be too small.
           </ul>

         @param nfev On exit, this is set to the number of function
           evaluations.

         @param wa Temporary storage array, of length <code>n</code>.
       */
      bool run(
        const lbfgs_parameters& lbfgs_params,
        int n,
        double* x,
        double f,
        double* g,
        double* s,
        int is0,
        double* stp,
        double ftol,
        double xtol,
        int maxfev,
        int* info,
        int* nfev,
        double* wa)
      {
        p5 = 0.5;
        p66 = 0.66;
        xtrapf = 4;
        if (info[0] != -1) {
          infoc[0] = 1;
          if (  n <= 0 || stp[0] <= 0 || ftol < 0
              || lbfgs_params.gtol < 0 || xtol < 0 || lbfgs_params.stpmin < 0
              || lbfgs_params.stpmax < lbfgs_params.stpmin || maxfev <= 0) {
            return true;
          }
          // Compute the initial gradient in the search direction
          // and check that s is a descent direction.
          dginit = 0;
          for (j = 1; j <= n; j += 1) {
            dginit = dginit + g [j -1] * s [is0+j -1];
          }
          if (dginit >= 0) {
            return false;
          }
          brackt[0] = false;
          stage1 = true;
          nfev[0] = 0;
          finit = f;
          dgtest = ftol*dginit;
          width = lbfgs_params.stpmax - lbfgs_params.stpmin;
          width1 = width/p5;
          for (j = 1; j <= n; j += 1) {
            wa [j -1] = x [j -1];
          }
          // The variables stx, fx, dgx contain the values of the step,
          // function, and directional derivative at the best step.
          // The variables sty, fy, dgy contain the value of the step,
          // function, and derivative at the other endpoint of
          // the interval of uncertainty.
          // The variables stp, f, dg contain the values of the step,
          // function, and derivative at the current step.
          stx[0] = 0;
          fx[0] = finit;
          dgx[0] = dginit;
          sty[0] = 0;
          fy[0] = finit;
          dgy[0] = dginit;
        }
        for (;;) {
          if (info[0] != -1) {
            // Set the minimum and maximum steps to correspond
            // to the present interval of uncertainty.
            if (brackt[0]) {
              stmin = std::min(stx[0], sty[0]);
              stmax = std::max(stx[0], sty[0]);
            }
            else {
              stmin = stx[0];
              stmax = stp[0] + xtrapf * (stp[0] - stx[0]);
            }
            // Force the step to be within the bounds stpmax and stpmin.
            stp[0] = std::max(stp[0], lbfgs_params.stpmin);
            stp[0] = std::min(stp[0], lbfgs_params.stpmax);
            // If an unusual termination is to occur then let
            // stp be the lowest point obtained so far.
            if (   (brackt[0] && (stp[0] <= stmin || stp[0] >= stmax))
                || nfev[0] >= maxfev - 1 || infoc[0] == 0
                || (brackt[0] && stmax - stmin <= xtol * stmax)) {
              stp[0] = stx[0];
            }
            // Evaluate the function and gradient at stp
            // and compute the directional derivative.
            // We return to main program to obtain F and G.
            for (j = 1; j <= n; j += 1) {
              x [j -1] = wa [j -1] + stp[0] * s [is0+j -1];
            }
            info[0]=-1;
            break;
          }
          info[0]=0;
          nfev[0] = nfev[0] + 1;
          dg = 0;
          for (j = 1; j <= n; j += 1) {
            dg = dg + g [j -1] * s [is0+j -1];
          }
          ftest1 = finit + stp[0]*dgtest;
          // Test for convergence.
          if (   (brackt[0] && (stp[0] <= stmin || stp[0] >= stmax))
              || infoc[0] == 0) {
            info[0] = 6;
          }
          if (stp[0] == lbfgs_params.stpmax && f <= ftest1 && dg <= dgtest) {
            info[0] = 5;
          }
          if (stp[0] == lbfgs_params.stpmin && (f > ftest1 || dg >= dgtest)) {
            info[0] = 4;
          }
          if (nfev[0] >= maxfev) {
            info[0] = 3;
          }
          if (brackt[0] && stmax - stmin <= xtol * stmax) {
            info[0] = 2;
          }
          if (f <= ftest1 && std::abs(dg) <= lbfgs_params.gtol * (-dginit)) {
            info[0] = 1;
          }
          // Check for termination.
          if (info[0] != 0) break;
          // In the first stage we seek a step for which the modified
          // function has a nonpositive value and nonnegative derivative.

          if (   stage1 && f <= ftest1
              && dg >= std::min(ftol, lbfgs_params.gtol) * dginit) {
            stage1 = false;
          }
          // A modified function is used to predict the step only if
          // we have not obtained a step for which the modified
          // function has a nonpositive function value and nonnegative
          // derivative, and if a lower function value has been
          // obtained but the decrease is not sufficient.
          if (stage1 && f <= fx[0] && f > ftest1) {
            // Define the modified function and derivative values.
            fm = f - stp[0]*dgtest;
            fxm[0] = fx[0] - stx[0]*dgtest;
            fym[0] = fy[0] - sty[0]*dgtest;
            dgm = dg - dgtest;
            dgxm[0] = dgx[0] - dgtest;
            dgym[0] = dgy[0] - dgtest;
            // Call cstep to update the interval of uncertainty
            // and to compute the new step.
            mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm,
                   brackt, stmin, stmax, infoc);
            // Reset the function and gradient values for f.
            fx[0] = fxm[0] + stx[0]*dgtest;
            fy[0] = fym[0] + sty[0]*dgtest;
            dgx[0] = dgxm[0] + dgtest;
            dgy[0] = dgym[0] + dgtest;
          }
          else {
            // Call mcstep to update the interval of uncertainty
            // and to compute the new step.
            mcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg,
                   brackt, stmin, stmax, infoc);
          }
          // Force a sufficient decrease in the size of the
          // interval of uncertainty.
          if (brackt[0]) {
            if (std::abs(sty[0] - stx[0]) >= p66 * width1) {
              stp[0] = stx[0] + p5 * (sty[0] - stx[0]);
            }
            width1 = width;
            width = std::abs(sty[0] - stx[0]);
          }
        }
        return true;
      }

      /* The purpose of this function is to compute a safeguarded step
         for a linesearch and to update an interval of uncertainty for
         a minimizer of the function.<p>

         The parameter <code>stx</code> contains the step with the
         least function value. The parameter <code>stp</code> contains
         the current step. It is assumed that the derivative at
         <code>stx</code> is negative in the direction of the step. If
         <code>brackt[0]</code> is <code>true</code> when
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
           been bracketed, then on output this variable is <code>true</code>.

         @param stpmin Lower bound for the step.
         @param stpmax Upper bound for the step.

         @param info On return from <code>mcstep</code>, this is set as
           follows:
           If <code>info</code> is 1, 2, 3, or 4, then the step has been
           computed successfully. Otherwise <code>info</code> = 0, and this
           indicates improper input parameters.

         @author Jorge J. More, David J. Thuente: original Fortran version,
           as part of Minpack project. Argonne Nat'l Laboratory, June 1983.
           Robert Dodier: Java translation, August 1997.
       */
      void mcstep(
        double* stx,
        double* fx,
        double* dx,
        double* sty,
        double* fy,
        double* dy,
        double* stp,
        double fp,
        double dp,
        bool* brackt,
        double stpmin,
        double stpmax,
        int* info)
      {
        bool bound;
        double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
        info[0] = 0;
        if (   (   brackt[0] && (stp[0] <= std::min(stx[0], sty[0])
                || stp[0] >= std::max(stx[0], sty[0])))
            || dx[0] * (stp[0] - stx[0]) >= 0.0 || stpmax < stpmin) {
          return;
        }
        // Determine if the derivatives have opposite sign.
        sgnd = dp * (dx[0] / std::abs(dx[0]));
        if (fp > fx[0]) {
          // First case. A higher function value.
          // The minimum is bracketed. If the cubic step is closer
          // to stx than the quadratic step, the cubic step is taken,
          // else the average of the cubic and quadratic steps is taken.
          info[0] = 1;
          bound = true;
          theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
          s = max3(std::abs(theta), std::abs(dx[0]), std::abs(dp));
          gamma = s * std::sqrt(sqr(theta / s) - (dx[0] / s) * (dp / s));
          if (stp[0] < stx[0]) gamma = - gamma;
          p = (gamma - dx[0]) + theta;
          q = ((gamma - dx[0]) + gamma) + dp;
          r = p/q;
          stpc = stx[0] + r * (stp[0] - stx[0]);
          stpq = stx[0]
            + ((dx[0] / ((fx[0] - fp) / (stp[0] - stx[0]) + dx[0])) / 2)
              * (stp[0] - stx[0]);
          if (std::abs(stpc - stx[0]) < std::abs(stpq - stx[0])) {
            stpf = stpc;
          }
          else {
            stpf = stpc + (stpq - stpc) / 2;
          }
          brackt[0] = true;
        }
        else if (sgnd < 0.0) {
          // Second case. A lower function value and derivatives of
          // opposite sign. The minimum is bracketed. If the cubic
          // step is closer to stx than the quadratic (secant) step,
          // the cubic step is taken, else the quadratic step is taken.
          info[0] = 2;
          bound = false;
          theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
          s = max3(std::abs(theta), std::abs(dx[0]), std::abs(dp));
          gamma = s * std::sqrt(sqr(theta / s) - (dx[0] / s) * (dp / s));
          if (stp[0] > stx[0]) gamma = - gamma;
          p = (gamma - dp) + theta;
          q = ((gamma - dp) + gamma) + dx[0];
          r = p/q;
          stpc = stp[0] + r * (stx[0] - stp[0]);
          stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
          if (std::abs(stpc - stp[0]) > std::abs(stpq - stp[0])) {
            stpf = stpc;
          }
          else {
            stpf = stpq;
          }
          brackt[0] = true;
        }
        else if (std::abs(dp) < std::abs(dx[0])) {
          // Third case. A lower function value, derivatives of the
          // same sign, and the magnitude of the derivative decreases.
          // The cubic step is only used if the cubic tends to infinity
          // in the direction of the step or if the minimum of the cubic
          // is beyond stp. Otherwise the cubic step is defined to be
          // either stpmin or stpmax. The quadratic (secant) step is also
          // computed and if the minimum is bracketed then the the step
          // closest to stx is taken, else the step farthest away is taken.
          info[0] = 3;
          bound = true;
          theta = 3 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
          s = max3(std::abs(theta), std::abs(dx[0]), std::abs(dp));
          gamma = s * std::sqrt(
            std::max(double(0), sqr(theta / s) - (dx[0] / s) * (dp / s)));
          if (stp[0] > stx[0]) gamma = -gamma;
          p = (gamma - dp) + theta;
          q = (gamma + (dx[0] - dp)) + gamma;
          r = p/q;
          if (r < 0.0 && gamma != 0.0) {
            stpc = stp[0] + r * (stx[0] - stp[0]);
          }
          else if (stp[0] > stx[0]) {
            stpc = stpmax;
          }
          else {
            stpc = stpmin;
          }
          stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
          if (brackt[0]) {
            if (std::abs(stp[0] - stpc) < std::abs(stp[0] - stpq)) {
              stpf = stpc;
            }
            else {
              stpf = stpq;
            }
          }
          else {
            if (std::abs(stp[0] - stpc) > std::abs(stp[0] - stpq)) {
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
          info[0] = 4;
          bound = false;
          if (brackt[0]) {
            theta = 3 * (fp - fy[0]) / (sty[0] - stp[0]) + dy[0] + dp;
            s = max3(std::abs(theta), std::abs(dy[0]), std::abs(dp));
            gamma = s * std::sqrt(sqr(theta / s) - (dy[0] / s) * (dp / s));
            if (stp[0] > sty[0]) gamma = -gamma;
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy[0];
            r = p/q;
            stpc = stp[0] + r * (sty[0] - stp[0]);
            stpf = stpc;
          }
          else if (stp[0] > stx[0]) {
            stpf = stpmax;
          }
          else {
            stpf = stpmin;
          }
        }
        // Update the interval of uncertainty. This update does not
        // depend on the new step or the case analysis above.
        if (fp > fx[0]) {
          sty[0] = stp[0];
          fy[0] = fp;
          dy[0] = dp;
        }
        else {
          if (sgnd < 0.0) {
            sty[0] = stx[0];
            fy[0] = fx[0];
            dy[0] = dx[0];
          }
          stx[0] = stp[0];
          fx[0] = fp;
          dx[0] = dp;
        }
        // Compute the new step and safeguard it.
        stpf = std::min(stpmax, stpf);
        stpf = std::max(stpmin, stpf);
        stp[0] = stpf;
        if (brackt[0] && bound) {
          if (sty[0] > stx[0]) {
            stp[0] = std::min(stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
          }
          else {
            stp[0] = std::max(stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
          }
        }
      }
  };

  class lbfgs
  {
    public:
      /*! <code>gtol</code> controls the accuracy of the line search
          <code>mcsrch</code>.
          If the function and gradient evaluations are inexpensive with
          respect to the cost of the iteration (which is sometimes the
          case when solving very large problems) it may be advantageous
          to set <code>gtol</code> to a small value. A typical small
          value is 0.1.
          Restriction: <code>gtol</code> should be greater than 1e-4.

          <code>stpmin</code> specifies the lower bound for the step
          in the line search.
          The default value is 1e-20. This value need not be modified
          unless the exponent is too large for the machine being used,
          or unless the problem is extremely badly scaled (in which
          case the exponent should be increased).

          <code>stpmax</code> specifies the upper bound for the step
          in the line search.
          The default value is 1e20. This value need not be modified
          unless the exponent is too large for the machine being used,
          or unless the problem is extremely badly scaled (in which
          case the exponent should be increased).
       */
      explicit lbfgs(
        double gtol = 0.9,
        double stpmin = 1.e-20,
        double stpmax = 1.e20)
        : params(gtol, stpmin, stpmax)
      {
        init();
      }

      //! Reinitialization.
      void init()
      {
        iflag = 0;
        gnorm = 0;
        stp1 = 0;
        ftol = 0;
        stp[0] = 0;
        ys = 0;
        yy = 0;
        sq = 0;
        yr = 0;
        beta = 0;
        xnorm = 0;
        iter = 0;
        nfun = 0;
        point = 0;
        ispt = 0;
        iypt = 0;
        maxfev = 0;
        info[0] = 0;
        bound = 0;
        npt = 0;
        cp = 0;
        i = 0;
        nfev[0] = 0;
        inmc = 0;
        iycn = 0;
        iscn = 0;
        finish = false;
      }

      /*! The solution vector as it was at the end of the most recently
          completed line search. This will usually be different from
          the return value of the parameter <tt>x</tt> of
          <tt>lbfgs</tt>, which is modified by line-search steps. A
          caller which wants to stop the optimization iterations before
          <tt>lbfgs</tt> automatically stops (by reaching a very
          small gradient) should copy this vector instead of using
          <tt>x</tt>. When <tt>lbfgs</tt> automatically stops,
          then <tt>x</tt> and <tt>solution_cache</tt> are the same.
        */
      std::vector<double> solution_cache;

      /*! This subroutine solves the unconstrained minimization problem
          <pre>
              min f(x),  x = (x1,x2,...,x_n),
          </pre>
          using the limited-memory BFGS method. The routine is
          especially effective on problems involving a large number of
          variables. In a typical iteration of this method an
          approximation <code>Hk</code> to the inverse of the Hessian
          is obtained by applying <code>m</code> BFGS updates to a
          diagonal matrix <code>Hk0</code>, using information from the
          previous M steps.  The user specifies the number
          <code>m</code>, which determines the amount of storage
          required by the routine. The user may also provide the
          diagonal matrices <code>Hk0</code> if not satisfied with the
          default choice.  The algorithm is described in "On the
          limited memory BFGS method for large scale optimization", by
          D. Liu and J. Nocedal, Mathematical Programming B 45 (1989)
          503-528.

          The user is required to calculate the function value
          <code>f</code> and its gradient <code>g</code>. In order to
          allow the user complete control over these computations,
          reverse  communication is used. The routine must be called
          repeatedly under the control of the parameter
          <code>iflag</code>.

          The steplength is determined at each iteration by means of
          the line search routine <code>mcsrch</code>, which is a
          slight modification of the routine <code>CSRCH</code> written
          by More' and Thuente.

          The only variables that are machine-dependent are
          <code>xtol</code>, <code>stpmin</code> and
          <code>stpmax</code>.

          Progress messages are accumulated in <code>log()</code>.
          Fatal errors cause exceptions to be thrown.

          @param n The number of variables in the minimization problem.
             Restriction: <code>n &gt; 0</code>.

          @param m The number of corrections used in the BFGS update.
             Values of <code>m</code> less than 3 are not recommended;
             large values of <code>m</code> will result in excessive
             computing time. <code>3 &lt;= m &lt;= 7</code> is
             recommended.
             Restriction: <code>m &gt; 0</code>.

          @param x On initial entry this must be set by the user to the
             values of the initial estimate of the solution vector. On
             exit with <code>iflag = 0</code>, it contains the values
             of the variables at the best point found (usually a
             solution).

          @param f Before initial entry and on a re-entry with
             <code>iflag = 1</code>, it must be set by the user to
             contain the value of the function <code>f</code> at the
             point <code>x</code>.

          @param g Before initial entry and on a re-entry with
             <code>iflag = 1</code>, it must be set by the user to
             contain the components of the gradient <code>g</code> at
             the point <code>x</code>.

          @param diagco Set this to <code>true</code> if the user
             wishes to provide the diagonal matrix <code>Hk0</code> at
             each iteration.  Otherwise it should be set to
             <code>false</code> in which case <code>lbfgs</code> will
             use a default value described below. If
             <code>diagco</code> is set to <code>true</code> the
             routine will return at each iteration of the algorithm
             with <code>iflag = 2</code>, and the diagonal matrix
             <code>Hk0</code> must be provided in the array
             <code>diag</code>.

          @param diag If <code>diagco = true</code>, then on initial
             entry or on re-entry with <code>iflag = 2</code>,
             <code>diag</code> must be set by the user to contain the
             values of the diagonal matrix <code>Hk0</code>.
             Restriction: all elements of <code>diag</code> must be
             positive.

          @param iprint Specifies output generated by <code>lbfgs</code>.
             <code>iprint[0]</code> specifies the frequency of the output:
             <ul>
             <li> <code>iprint[0] &lt; 0</code>: no output is generated,
             <li> <code>iprint[0] = 0</code>: output only at first and
                  last iteration,
             <li> <code>iprint[0] &gt; 0</code>: output every
                  <code>iprint[0]</code> iterations.
             </ul>

             <code>iprint[1]</code> specifies the type of output generated:
             <ul>
             <li><code>iprint[1] = 0</code>: iteration count, number of
                 function evaluations, function value, norm of the gradient,
                 and steplength,
             <li><code>iprint[1] = 1</code>: same as <code>iprint[1]=0</code>,
                 plus vector of variables and  gradient vector at the initial
                 point,
             <li><code>iprint[1] = 2</code>: same as <code>iprint[1]=1</code>,
                 plus vector of variables,
             <li><code>iprint[1] = 3</code>: same as <code>iprint[1]=2</code>,
                 plus gradient vector.
             </ul>

           @param eps Determines the accuracy with which the solution
             is to be found. The subroutine terminates when
             <pre>
                     ||G|| &lt; eps max(1,||X||),
             </pre>
             where <code>||.||</code> denotes the Euclidean norm.


           @param xtol An estimate of the machine precision (e.g. 10e-16
             on a SUN station 3/60). The line search routine will
             terminate if the relative width of the interval of
             uncertainty is less than <code>xtol</code>.

           @param iflag This must be set to 0 on initial entry to
             <code>lbfgs</code>.
             A return with <code>iflag &lt; 0</code> indicates an error,
             and <code>iflag = 0</code> indicates that the routine has
             terminated without detecting errors. On a return with
             <code>iflag = 1</code>, the user must evaluate the function
             <code>f</code> and gradient <code>g</code>. On a return with
             <code>iflag = 2</code>, the user must provide the diagonal
             matrix <code>Hk0</code>.

             The following negative values of <code>iflag</code>,
             detecting an error, are possible:

             <ul>
             <li><code>iflag = -1</code> The line search routine
                 <code>mcsrch</code> failed. One of the following
                 messages is printed:
               <ul>
               <li>Improper input parameters.
               <li>Relative width of the interval of uncertainty is at
                   most <code>xtol</code>.
               <li>More than 20 function evaluations were required at the
                   present iteration.
               <li>The step is too small.
               <li>The step is too large.
               <li>Rounding errors prevent further progress. There may not
                   be a step which satisfies the sufficient decrease and
                   curvature conditions. Tolerances may be too small.
               </ul>
             <li><code>iflag = -2</code> The i-th diagonal element of
                 the diagonal inverse Hessian approximation, given in DIAG,
                 is not positive.
             <li><code>iflag = -3</code> Improper input parameters for LBFGS
                 (<code>n</code> or <code>m</code> are not positive).
             </ul>
       */
      void run(
        int n,
        int m,
        double* x,
        double f,
        double* g,
        bool diagco,
        double* diag,
        const int* iprint,
        double eps,
        double xtol)
      {
        bool execute_entire_while_loop = false;
        m_w.resize(n*(2*m+1)+2*m);
        double* w = &(*(m_w.begin()));
        if (iflag == 0) { // Initialize.
          solution_cache.assign(x, x + n);
          iter = 0;
          if (n <= 0 || m <= 0) {
            iflag = -3;
            throw error(
              "Improper input parameters (n or m are not positive.)");
          }
          if (params.gtol <= 0.0001) {
            m_log.push_back(
              "gtol is less than or equal to 0.0001."
              " It has been reset to 0.9.");
            params.gtol= 0.9;
          }
          nfun= 1;
          point= 0;
          finish= false;
          if (diagco) {
            for (i = 1; i <= n; i += 1) {
              if (diag[i -1] <= 0) {
                iflag = -2;
                throw_diagonal_element_not_positive(i);
              }
            }
          }
          else {
            for (i = 1; i <= n; i += 1) {
              diag[i -1] = 1;
            }
          }
          ispt= n+2*m;
          iypt= ispt+n*m;
          for (i = 1; i <= n; i += 1) {
            w[ispt + i -1] = -g[i -1] * diag[i -1];
          }
          gnorm = std::sqrt(ddot(n, g, 0, 1, g, 0, 1));
          stp1= 1/gnorm;
          ftol= 0.0001;
          maxfev= 20;
          if (iprint[1 -1] >= 0) {
            lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
          }
          execute_entire_while_loop = true;
        }
        for (;;) {
          if (execute_entire_while_loop) {
            iter= iter+1;
            info[0]=0;
            bound=iter-1;
            if (iter != 1) {
              if (iter > m) bound = m;
              ys = ddot(n, w, iypt + npt, 1, w, ispt + npt, 1);
              if (!diagco) {
                yy = ddot(n, w, iypt + npt, 1, w, iypt + npt, 1);
                for (i = 1; i <= n; i += 1) {
                  diag[i -1] = ys / yy;
                }
              }
              else {
                iflag = 2;
                return;
              }
            }
          }
          if (execute_entire_while_loop || iflag == 2) {
            if (iter != 1) {
              if (diagco) {
                for (i = 1; i <= n; i += 1) {
                  if (diag[i -1] <= 0) {
                    iflag = -2;
                    throw_diagonal_element_not_positive(i);
                  }
                }
              }
              cp = point;
              if (point == 0) cp = m;
              w[n + cp -1] = 1 / ys;
              for (i = 1; i <= n; i += 1) {
                w[i -1] = -g[i -1];
              }
              cp = point;
              for (i = 1; i <= bound; i += 1) {
                cp=cp-1;
                if (cp == -1) cp = m - 1;
                sq = ddot(n, w, ispt + cp * n, 1, w, 0, 1);
                inmc=n+m+cp+1;
                iycn=iypt+cp*n;
                w[inmc -1] = w[n + cp + 1 -1] * sq;
                daxpy(n, -w[inmc -1], w, iycn, 1, w, 0, 1);
              }
              for (i = 1; i <= n; i += 1) {
                w[i -1] = diag[i -1] * w[i -1];
              }
              for (i = 1; i <= bound; i += 1) {
                yr = ddot(n, w, iypt + cp * n, 1, w, 0, 1);
                beta = w[n + cp + 1 -1] * yr;
                inmc=n+m+cp+1;
                beta = w[inmc -1] - beta;
                iscn=ispt+cp*n;
                daxpy(n, beta, w, iscn, 1, w, 0, 1);
                cp=cp+1;
                if (cp == m) cp = 0;
              }
              for (i = 1; i <= n; i += 1) {
                w[ispt + point * n + i -1] = w[i -1];
              }
            }
            nfev[0]=0;
            stp[0]=1;
            if (iter == 1) stp[0] = stp1;
            for (i = 1; i <= n; i += 1) {
              w[i -1] = g[i -1];
            }
          }
          if (!mcsrch_instance.run(params, n, x, f, g, w, ispt + point * n,
                stp, ftol, xtol, maxfev, info, nfev, diag)) {
            m_log.push_back(
              "The search direction is not a descent direction.\n");
          }
          if (info[0] == -1) {
            iflag = 1;
            return;
          }
          if (info[0] != 1) {
            iflag = -1;
            throw error(
              "Line search failed. See documentation of routine mcsrch."
              " Error return of line search: info = " + itoa(info[0]) +
              " Possible causes: function or gradient are incorrect,"
              " or incorrect tolerances.");
          }
          nfun= nfun + nfev[0];
          npt=point*n;
          for (i = 1; i <= n; i += 1) {
            w[ispt + npt + i -1] = stp[0] * w[ispt + npt + i -1];
            w[iypt + npt + i -1] = g[i -1] - w[i -1];
          }
          point=point+1;
          if (point == m) point = 0;
          gnorm = std::sqrt(ddot(n, g, 0, 1, g, 0, 1));
          xnorm = std::sqrt(ddot(n, x, 0, 1, x, 0, 1));
          xnorm = std::max(1.0, xnorm);
          if (gnorm / xnorm <= eps) finish = true;
          if (iprint[1 -1] >= 0) {
            lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
          }
          // Cache the current solution vector. Due to the spaghetti-like
          // nature of this code, it's not possible to quit here and return;
          // we need to go back to the top of the loop, and eventually call
          // mcsrch one more time -- but that will modify the solution vector.
          // So we need to keep a copy of the solution vector as it was at
          // the completion (info[0]==1) of the most recent line search.
          solution_cache.assign(x, x + n);
          if (finish) {
            iflag = 0;
            return;
          }
          execute_entire_while_loop = true; // from now on, execute whole loop
        }
      }

      /*! This method returns the total number of evaluations of the
          objective function since the last time LBFGS was restarted.
          The total number of function evaluations increases by the
          number of evaluations required for the line search; the total
          is only increased after a successful line search.
        */
      int nfevaluations() const { return nfun; }

      int flag() const { return iflag; }

      const af::shared<std::string>& log() const { return m_log; }
            af::shared<std::string>  log()       { return m_log; }

    protected:

      static void throw_diagonal_element_not_positive(int i) {
        throw error(
          "The " + itoa(i) + "-th diagonal element of the inverse"
          " Hessian approximation is not positive.");
      }

      /*! Print debugging and status messages for <code>lbfgs</code>.
          Depending on the parameter <code>iprint</code>, this can include
          number of function evaluations, current function value, etc.
          The messages are accumulated in <code>log()</code>.

          @param iprint Specifies output generated by <code>lbfgs</code>.<p>
             <code>iprint[0]</code> specifies the frequency of the output:
             <ul>
             <li><code>iprint[0] &lt; 0</code>:
                 no output is generated,
             <li><code>iprint[0] = 0</code>:
                 output only at first and last iteration,
             <li><code>iprint[0] &gt; 0</code>:
                 output every <code>iprint[0]</code> iterations.
             </ul>
             <p>
             <code>iprint[1]</code> specifies the type of output generated:
             <ul>
             <li><code>iprint[1] = 0</code>:
                 iteration count, number of function evaluations,
                 function value, norm of the gradient, and steplength,
             <li><code>iprint[1] = 1</code>:
                 same as <code>iprint[1]=0</code>, plus vector of
                 variables and  gradient vector at the initial point,
             <li><code>iprint[1] = 2</code>:
                 same as <code>iprint[1]=1</code>, plus vector of variables,
             <li><code>iprint[1] = 3</code>:
                 same as <code>iprint[1]=2</code>, plus gradient vector.
             </ul>
          @param iter Number of iterations so far.
          @param nfun Number of function evaluations so far.
          @param gnorm Norm of gradient at current solution <code>x</code>.
          @param n Number of free parameters.
          @param m Number of corrections kept.
          @param x Current solution.
          @param f Function value at current solution.
          @param g Gradient at current solution <code>x</code>.
          @param stp Current stepsize.
          @param finish Whether this method should print
                 the ``we're done'' message.
        */
      void lb1(
        const int* iprint,
        int iter,
        int nfun,
        double gnorm,
        int n,
        int m,
        double* x,
        double f,
        double* g,
        double* stp,
        bool finish)
      {
        static const char* header =
        "        i       nfn           func           gnorm      steplength\n";
        int i;
        if (iter == 0) {
          m_log.push_back(
            "*************************************************\n");
          m_log.push_back(
            "  n = " + itoa(n) + "   number of corrections = "
            + m + "\n       initial values\n");
          m_log.push_back(
            " f =  " + ftoa(f) + "   gnorm =  " + gnorm + "\n");
          if (iprint[2 -1] >= 1) {
            m_log.push_back(" vector x =");
            for (i = 1; i <= n; i++) {
              m_log.push_back("  " + ftoa(x[i-1]));
            }
            m_log.push_back("\n");
            m_log.push_back(" gradient vector g =");
            for (i = 1; i <= n; i++) {
              m_log.push_back("  " + ftoa(g[i-1]));
            }
            m_log.push_back("\n");
          }
          m_log.push_back(
            "*************************************************\n");
          m_log.push_back(header);
        }
        else {
          if ((iprint[1 -1] == 0) && (iter != 1 && ! finish)) return;
          if (iprint[1 -1] != 0) {
            if ((iter - 1) % iprint[1 -1] == 0 || finish) {
              if (iprint[2 -1] > 1 && iter > 1) {
                m_log.push_back(header);
              }
              m_log.push_back(
                "\t"+itoa(iter)+"\t"+nfun+"\t"+f+"\t"+gnorm+"\t"+stp[0] + "\n");
            }
            else {
              return;
            }
          }
          else {
            if (iprint[2 -1] > 1 && finish) {
              m_log.push_back(header);
            }
            m_log.push_back(
              "\t"+itoa(iter)+"\t"+nfun+"\t"+f+"\t"+gnorm+"\t"+stp[0] + "\n");
          }
          if (iprint[2 -1] == 2 || iprint[2 -1] == 3) {
            if (finish) {
              m_log.push_back(" final point x =");
            }
            else {
              m_log.push_back(" vector x =  ");
            }
            for (i = 1; i <= n; i++) {
              m_log.push_back("  " + ftoa(x[i-1]));
            }
            m_log.push_back("\n");
            if (iprint[2 -1] == 3) {
              m_log.push_back(" gradient vector g =");
              for (i = 1; i <= n; i++) {
                m_log.push_back("  " + ftoa(g[i-1]));
              }
              m_log.push_back("\n");
            }
          }
          if (finish) {
            m_log.push_back(
              " The minimization terminated without detecting errors.\n"
              " iflag = 0\n");
          }
        }
      }

      /*! Compute the sum of a vector times a scalara plus another vector.
          Adapted from the subroutine <code>daxpy</code> in
          <code>lbfgs.f</code>.
          There could well be faster ways to carry out this operation; this
          code is a straight translation from the Fortran.
        */
      void daxpy(
        int n,
        double da,
        double* dx,
        int ix0,
        int incx,
        double* dy,
        int iy0,
        int incy)
      {
        int i, ix, iy, m, mp1;
        if (n <= 0) return;
        if (da == 0) return;
        if  (!(incx == 1 && incy == 1)) {
          ix = 1;
          iy = 1;
          if (incx < 0) ix = (-n + 1) * incx + 1;
          if (incy < 0) iy = (-n + 1) * incy + 1;
          for (i = 1; i <= n; i += 1) {
            dy[iy0+iy -1] = dy[iy0+iy -1] + da * dx[ix0+ix -1];
            ix = ix + incx;
            iy = iy + incy;
          }
          return;
        }
        m = n % 4;
        if (m != 0) {
          for (i = 1; i <= m; i += 1) {
            dy[iy0+i -1] = dy[iy0+i -1] + da * dx[ix0+i -1];
          }
          if (n < 4) return;
        }
        mp1 = m + 1;
        for (i = mp1; i <= n; i += 4) {
          dy[iy0+i -1] = dy[iy0+i -1] + da * dx[ix0+i -1];
          dy[iy0+i + 1 -1] = dy[iy0+i + 1 -1] + da * dx[ix0+i + 1 -1];
          dy[iy0+i + 2 -1] = dy[iy0+i + 2 -1] + da * dx[ix0+i + 2 -1];
          dy[iy0+i + 3 -1] = dy[iy0+i + 3 -1] + da * dx[ix0+i + 3 -1];
        }
      }

      /*! Compute the dot product of two vectors.
          Adapted from the subroutine <code>ddot</code>
          in <code>lbfgs.f</code>.
          There could well be faster ways to carry out this operation;
          this code is a straight translation from the Fortran.
        */
      double ddot(
        int n,
        double* dx,
        int ix0,
        int incx,
        double* dy,
        int iy0,
        int incy)
      {
        int i, ix, iy, m, mp1;
        double dtemp(0);
        if (n <= 0) return 0;
        if (!(incx == 1 && incy == 1)) {
          ix = 1;
          iy = 1;
          if (incx < 0) ix = (-n + 1) * incx + 1;
          if (incy < 0) iy = (-n + 1) * incy + 1;
          for (i = 1; i <= n; i += 1) {
            dtemp = dtemp + dx[ix0+ix -1] * dy[iy0+iy -1];
            ix = ix + incx;
            iy = iy + incy;
          }
          return dtemp;
        }
        m = n % 5;
        if (m != 0) {
          for (i = 1; i <= m; i += 1) {
            dtemp = dtemp + dx[ix0+i -1] * dy[iy0+i -1];
          }
          if (n < 5) return dtemp;
        }
        mp1 = m + 1;
        for (i = mp1; i <= n; i += 5) {
          dtemp += dx[ix0+i -1] * dy[iy0+i -1]
                 + dx[ix0+i + 1 -1] * dy[iy0+i + 1 -1]
                 + dx[ix0+i + 2 -1] * dy[iy0+i + 2 -1]
                 + dx[ix0+i + 3 -1] * dy[iy0+i + 3 -1]
                 + dx[ix0+i + 4 -1] * dy[iy0+i + 4 -1];
        }
        return dtemp;
      }
      lbfgs_parameters params;
      mcsrch mcsrch_instance;
      int iflag;
      double gnorm;
      double stp1;
      double ftol;
      double stp[1];
      double ys;
      double yy;
      double sq;
      double yr;
      double beta;
      double xnorm;
      int iter;
      int nfun;
      int point;
      int ispt;
      int iypt;
      int maxfev;
      int info[1];
      int bound;
      int npt;
      int cp;
      int i;
      int nfev[1];
      int inmc;
      int iycn;
      int iscn;
      bool finish;
      std::vector<double> m_w;
      af::shared<std::string> m_log;
  };

} // namespace cctbx

#endif // CCTBX_LBFGS_H
