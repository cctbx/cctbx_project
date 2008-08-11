#ifndef SCITBX_LINE_SEARCH_MORE_THUENTE_1994_RAW_H
#define SCITBX_LINE_SEARCH_MORE_THUENTE_1994_RAW_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <cstddef>

namespace scitbx { namespace line_search {

  template <typename FloatType=double>
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
      std::vector<FloatType> initial_x;

      static
      FloatType
      pow2(FloatType const& x) { return x * x; }

      static
      FloatType const&
      max3(
        FloatType const& x,
        FloatType const& y,
        FloatType const& z)
      {
        return x < y ? (y < z ? z : y ) : (x < z ? z : x );
      }

    public:
      const char* info_meaning;

      void
      free_workspace() { initial_x = std::vector<FloatType>(); }

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
       */
      void
      run(
        FloatType const& gtol,
        FloatType const& stpmin,
        FloatType const& stpmax,
        unsigned n,
        FloatType* x,
        FloatType f,
        const FloatType* g,
        const FloatType* s,
        FloatType& stp,
        FloatType ftol,
        FloatType xtol,
        unsigned maxfev,
        int& info,
        unsigned& nfev);

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
      static
      int
      mcstep(
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
    FloatType const& gtol,
    FloatType const& stpmin,
    FloatType const& stpmax,
    unsigned n,
    FloatType* x,
    FloatType f,
    const FloatType* g,
    const FloatType* s,
    FloatType& stp,
    FloatType ftol,
    FloatType xtol,
    unsigned maxfev,
    int& info,
    unsigned& nfev)
  {
    if (info != -1) {
      infoc = 1;
      if (   n == 0
          || maxfev == 0
          || gtol < FloatType(0)
          || xtol < FloatType(0)
          || stpmin < FloatType(0)
          || stpmax < stpmin) {
        throw std::runtime_error("Improper input parameters.");
      }
      if (stp <= FloatType(0) || ftol < FloatType(0)) {
        throw std::runtime_error("Improper value for stp or ftol.");
      }
      // Compute the initial gradient in the search direction
      // and check that s is a descent direction.
      dginit = FloatType(0);
      for (unsigned j = 0; j < n; j++) {
        dginit += g[j] * s[j];
      }
      if (dginit >= FloatType(0)) {
        throw std::runtime_error("Search direction not descent.");
      }
      brackt = false;
      stage1 = true;
      nfev = 0;
      finit = f;
      dgtest = ftol*dginit;
      width = stpmax - stpmin;
      width1 = FloatType(2) * width;
      initial_x.assign(x, x+n);
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
        for (unsigned j = 0; j < n; j++) {
          x[j] = initial_x[j] + stp * s[j];
        }
        info = -1;
        info_meaning =
          "A return is made to compute the function and gradient.";
        break;
      }
      info = 0;
      info_meaning = 0;
      nfev++;
      FloatType dg(0);
      for (unsigned j = 0; j < n; j++) {
        dg += g[j] * s[j];
      }
      FloatType ftest1 = finit + stp*dgtest;
      // Test for convergence.
      if ((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0) {
        info = 6;
        info_meaning =
          "Rounding errors prevent further progress."
          " There may not be a step which satisfies the"
          " sufficient decrease and curvature conditions."
          " Tolerances may be too small.";
        break;
      }
      if (stp == stpmax && f <= ftest1 && dg <= dgtest) {
        info = 5;
        info_meaning = "The step is at the upper bound stpmax.";
        break;
      }
      if (stp == stpmin && (f > ftest1 || dg >= dgtest)) {
        info = 4;
        info_meaning = "The step is at the lower bound stpmin.";
        break;
      }
      if (nfev >= maxfev) {
        info = 3;
        info_meaning = "Number of function evaluations has reached maxfev.";
        break;
      }
      if (brackt && stmax - stmin <= xtol * stmax) {
        info = 2;
        info_meaning =
          "Relative width of the interval of uncertainty"
          " is at most xtol.";
        break;
      }
      // Check for termination.
      if (f <= ftest1 && std::abs(dg) <= gtol * (-dginit)) {
        info = 1;
        info_meaning =
          "The sufficient decrease condition and the"
          " directional derivative condition hold.";
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
        if (std::abs(sty - stx) >= FloatType(0.66) * width1) {
          stp = stx + FloatType(0.5) * (sty - stx);
        }
        width1 = width;
        width = std::abs(sty - stx);
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
    sgnd = dp * (dx / std::abs(dx));
    if (fp > fx) {
      // First case. A higher function value.
      // The minimum is bracketed. If the cubic step is closer
      // to stx than the quadratic step, the cubic step is taken,
      // else the average of the cubic and quadratic steps is taken.
      info = 1;
      bound = true;
      theta = FloatType(3) * (fx - fp) / (stp - stx) + dx + dp;
      s = max3(std::abs(theta), std::abs(dx), std::abs(dp));
      gamma = s * std::sqrt(pow2(theta / s) - (dx / s) * (dp / s));
      if (stp < stx) gamma = - gamma;
      p = (gamma - dx) + theta;
      q = ((gamma - dx) + gamma) + dp;
      r = p/q;
      stpc = stx + r * (stp - stx);
      stpq = stx
        + ((dx / ((fx - fp) / (stp - stx) + dx)) / FloatType(2))
          * (stp - stx);
      if (std::abs(stpc - stx) < std::abs(stpq - stx)) {
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
      s = max3(std::abs(theta), std::abs(dx), std::abs(dp));
      gamma = s * std::sqrt(pow2(theta / s) - (dx / s) * (dp / s));
      if (stp > stx) gamma = - gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dx;
      r = p/q;
      stpc = stp + r * (stx - stp);
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (std::abs(stpc - stp) > std::abs(stpq - stp)) {
        stpf = stpc;
      }
      else {
        stpf = stpq;
      }
      brackt = true;
    }
    else if (std::abs(dp) < std::abs(dx)) {
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
      s = max3(std::abs(theta), std::abs(dx), std::abs(dp));
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
        if (std::abs(stp - stpc) < std::abs(stp - stpq)) {
          stpf = stpc;
        }
        else {
          stpf = stpq;
        }
      }
      else {
        if (std::abs(stp - stpc) > std::abs(stp - stpq)) {
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
        s = max3(std::abs(theta), std::abs(dy), std::abs(dp));
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

}} // namespace scitbx::line_search

#endif // SCITBX_LINE_SEARCH_MORE_THUENTE_1994_RAW_H
