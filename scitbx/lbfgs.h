#ifndef SCITBX_LBFGS_H
#define SCITBX_LBFGS_H

#include <stdio.h>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <string>

namespace scitbx {

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
      error(std::string const& msg) throw()
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
        snprintf(buf, sizeof(buf), "%lu", i); // FUTURE: use C++ facility
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
      error_improper_input_parameter(std::string const& msg) throw()
        : error("Improper input parameter: " + msg)
      {}
  };

  //! Specific exception class.
  class error_improper_input_data : public error {
    public:
      //! Constructor.
      error_improper_input_data(std::string const& msg) throw()
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
      error_line_search_failed(std::string const& msg) throw()
        : error("Line search failed: " + msg)
      {}
  };

  //! Specific exception class.
  class error_line_search_failed_rounding_errors
  : public error_line_search_failed {
    public:
      //! Constructor.
      error_line_search_failed_rounding_errors(std::string const& msg) throw()
        : error_line_search_failed(msg)
      {}
  };

  namespace detail {

    /* Compute the sum of a vector times a scalar plus another vector.
       Adapted from the subroutine <code>daxpy</code> in
       <code>lbfgs.f</code>.
     */
    template <typename FloatType, typename SizeType>
    void daxpy(
      SizeType n,
      FloatType da,
      const FloatType* dx,
      SizeType ix0,
      SizeType incx,
      FloatType* dy,
      SizeType iy0,
      SizeType incy)
    {
      SizeType i, ix, iy, m;
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

    template <typename FloatType, typename SizeType>
    inline
    void daxpy(
      SizeType n,
      FloatType da,
      const FloatType* dx,
      SizeType ix0,
      FloatType* dy)
    {
      daxpy(n, da, dx, ix0, SizeType(1), dy, SizeType(0), SizeType(1));
    }

    /* Compute the dot product of two vectors.
       Adapted from the subroutine <code>ddot</code>
       in <code>lbfgs.f</code>.
     */
    template <typename FloatType, typename SizeType>
    FloatType ddot(
      SizeType n,
      const FloatType* dx,
      SizeType ix0,
      SizeType incx,
      const FloatType* dy,
      SizeType iy0,
      SizeType incy)
    {
      SizeType i, ix, iy, m;
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

    template <typename FloatType, typename SizeType>
    inline
    FloatType ddot(
      SizeType n,
      const FloatType* dx,
      const FloatType* dy)
    {
      return ddot(
        n, dx, SizeType(0), SizeType(1), dy, SizeType(0), SizeType(1));
    }

    template <typename NumType>
    inline
    NumType
    pow2(NumType const& x) { return x * x; }

    template <typename NumType>
    inline
    NumType
    abs(NumType const& x) {
      if (x < NumType(0)) return -x;
      return x;
    }

    // This class implements an algorithm for multi-dimensional line search.
    template <typename FloatType, typename SizeType = std::size_t>
    class mcsrch
    {
      protected:
        int infoc;
        bool gradient_only;
        bool line_search;
        FloatType dginit;
        bool brackt;
        bool stage1;
        FloatType finit;
        FloatType dgtest;
        FloatType width;
        FloatType width1;
        FloatType stx;
        std::vector<FloatType> sx;
        FloatType fx;
        FloatType dgx;
        std::vector<FloatType> gx;
        FloatType sty;
        FloatType fy;
        FloatType dgy;
        std::vector<FloatType> gy;
        FloatType stmin;
        FloatType stmax;

        static FloatType const& max3(
          FloatType const& x,
          FloatType const& y,
          FloatType const& z)
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
          bool gradient_only,
          bool line_search,
          FloatType const& gtol,
          FloatType const& stpmin,
          FloatType const& stpmax,
          SizeType n,
          FloatType* x,
          FloatType f,
          const FloatType* g,
          FloatType* s,
          SizeType is0,
          FloatType& stp,
          FloatType ftol,
          FloatType xtol,
          SizeType maxfev,
          int& info,
          SizeType& nfev,
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

        //---> Insertion starts
        static int mcstep_g(
         SizeType& n,
         FloatType* const& sx,
         FloatType* const& gx,
         FloatType& stp,
         FloatType* const& gp,
         bool& brackt,
         FloatType stpmin,
         FloatType stpmax);
         //<--- Insertion ends
    };

    template <typename FloatType, typename SizeType>
    void mcsrch<FloatType, SizeType>::run(
      bool gradient_only,
      bool line_search,
      FloatType const& gtol,
      FloatType const& stpmin,
      FloatType const& stpmax,
      SizeType n,
      FloatType* x,
      FloatType f,
      const FloatType* g,
      FloatType* s,
      SizeType is0,
      FloatType& stp,
      FloatType ftol,
      FloatType xtol,
      SizeType maxfev,
      int& info,
      SizeType& nfev,
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
        for (SizeType j = 0; j < n; j++) {
          dginit += g[j] * s[is0+j];
        }
        if (dginit >= FloatType(0) && !gradient_only) {
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
        gx = std::vector<FloatType>(n);
        sx = std::vector<FloatType>(n);
        for (SizeType j = 0; j < n; j++) {
          gx[j] = g[j];
          sx[j] = s[j+is0];
        }
        sty = FloatType(0);
        fy = finit;
        dgy = dginit;
      }
      for (;;) {
        if (info != -1) {
          //---> Insertion starts
          if (gradient_only){
            FloatType dr_max=0.0;
            SizeType m = n/3;
            for (SizeType i=0; i<m; i++){
              FloatType dr=0.0;
              SizeType j = 3*i;
              dr = sx[j]*sx[j];
              dr = dr + sx[j+1]*sx[j+1];
              dr = dr + sx[j+2]*sx[j+2];
              if(dr_max<dr){
                dr_max = dr;
              }
            }
            dr_max = std::sqrt(dr_max);
            if (abs(dr_max*stp) > stpmax) {
              stp = stpmax/dr_max;
            }
            else if (abs(dr_max*stp) < FloatType(0.0001)) {
              stp = FloatType(0.0001)/dr_max;
            }
          }
          //<--- Insertion ends
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
          for (SizeType j = 0; j < n; j++) {
            x[j] = wa[j] + stp * s[is0+j];
          }
          info=-1;
          break;
        }
        info = 0;
        nfev++;
        FloatType dg(0);
        for (SizeType j = 0; j < n; j++) {
          dg += g[j] * s[is0+j];
        }
        FloatType ftest1 = finit + stp*dgtest;
        // Check for termination.
        //---> Insertion starts
        if (gradient_only) {
          if (!line_search || nfev==2 || ((abs(dg) <= gtol * (-dginit)))) {
            info = 1;
            break;
          }
        }
        //<--- Insertion ends
        else {
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
            if  ( f <= ftest1 && abs(dg) <= gtol * (-dginit)) {
            info = 1;
            break;
            }
        }
        // In the first stage we seek a step for which the modified
        // function has a nonpositive value and nonnegative derivative.
        if (   stage1 && (f <= ftest1 || gradient_only)
            && dg >= std::min(ftol, gtol) * dginit) {
          stage1 = false;
        }
        //---> Insertion starts
        if (gradient_only) {
          gy = std::vector<FloatType>(n);
          for (SizeType j = 0; j < n; j++) {
            gy[j] = g[j];
          }
          FloatType* const& const_gy = &(gy[0]);
          FloatType* const& const_gx = &(gx[0]);
          FloatType* const& const_sx = &(sx[0]);
          infoc = mcstep_g(n,  const_sx, const_gx, stp,  const_gy,
                        brackt, stmin, stpmax);
         }
         //<--- Insertion ends
        else{
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
      } // end of for loop in mcsrch
    } // end of mcsrch

    // energy_gradient line search
    template <typename FloatType, typename SizeType>
    int mcsrch<FloatType, SizeType>::mcstep(
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

    // gradient_only line search
    //---> Insertion starts
    template <typename FloatType, typename SizeType>
    int mcsrch<FloatType, SizeType>::mcstep_g(
      SizeType& n,
      FloatType* const& sx,
      FloatType* const& gx,
      FloatType& stp,
      FloatType* const& gp,
      bool& brackt,
      FloatType stpmin,
      FloatType stpmax)
    {
      int info = 5;
      brackt = false;
      FloatType sxnorm = std::sqrt(ddot(SizeType(n), sx, sx));
      std::vector<FloatType> pk(n);
      std::vector<FloatType> a_sum(n);
      std::vector<FloatType> akm1(n);
      for (SizeType i=0; i < n; i++) {
        pk[i] = sx[i]/sxnorm;
        akm1[i] = -gx[i]/sxnorm;
        a_sum[i] = akm1[i] - gp[i]/sxnorm;
      }
      FloatType* const& a_sum_c = &(a_sum[0]);
      FloatType* const& akm1_c = &(akm1[0]);
      FloatType* const& pk_c = &(pk[0]);
      FloatType dfk = -0.5*ddot(SizeType(n), pk_c, a_sum_c);
      FloatType rho = -1.0*ddot(SizeType(n), pk_c, akm1_c);
      if (abs(dfk - rho) < 1e-5){
       rho = -(dfk - rho);
      }
      FloatType theta = -rho / (dfk - rho);
      if ((dfk < 0) && (theta < 0)){
        theta = -theta;
      }
      stp = theta / 2.0;
      return info;
    }
    //<--- Insertion ends
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
  template <typename FloatType, typename SizeType = std::size_t>
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

          @param maxfev Maximum number of function evaluations
             <b>per line search</b>.
             Termination occurs when the number of evaluations
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
        SizeType n,
        SizeType m = 5,
        SizeType maxfev = 20,
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
        scratch_array_.resize(n_);
      }

      //! Number of free parameters (as passed to the constructor).
      SizeType n() const { return n_; }

      //! Number of corrections kept (as passed to the constructor).
      SizeType m() const { return m_; }

      /*! \brief Maximum number of evaluations of the objective function
          per line search (as passed to the constructor).
       */
      SizeType maxfev() const { return maxfev_; }

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
      SizeType iter() const { return iter_; }

      //! Total number of evaluations of the objective function so far.
      /*! The total number of function evaluations increases by the
          number of evaluations required for the line search. The total
          is only increased after a successful line search.
          <p>
          See also: iter()
       */
      SizeType nfun() const { return nfun_; }

      //! Norm of gradient given gradient array of length n().
      FloatType euclidean_norm(const FloatType* a) const {
        return std::sqrt(detail::ddot(n_, a, a));
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
        return generic_run(x, f, g, false, false, true, 0);
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
        const FloatType* diag)
      //---> Insertion starts
      {
        return generic_run(x, f, g, true, false, true, diag);
      }
      //<--- Insertion ends

      bool run(
        FloatType* x,
        FloatType f,
        const FloatType* g,
      //---> Insertion starts
        bool gradient_only,
        bool line_search)
      {
        return generic_run(x, f, g, false, gradient_only, line_search, 0);
      }
      //<--- Insertion ends

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
        const FloatType* diag,
        //---> Insertion starts
        bool gradient_only,
        bool line_search)
      {
        return generic_run(x, f, g, true, gradient_only, line_search, diag);
      }
      //<--- Insertion ends

    protected:
      static void throw_diagonal_element_not_positive(SizeType i) {
        throw error_improper_input_data(
          "The " + error::itoa(i) + ". diagonal element of the"
          " inverse Hessian approximation is not positive.");
      }

      bool generic_run(
        FloatType* x,
        FloatType f,
        const FloatType* g,
        bool diagco,
        //---> Insertion starts
        bool gradient_only,
        bool line_search,
        //<--- Insertion ends
        const FloatType* diag);

      detail::mcsrch<FloatType, SizeType> mcsrch_instance;
      const SizeType n_;
      const SizeType m_;
      const SizeType maxfev_;
      const FloatType gtol_;
      const FloatType xtol_;
      const FloatType stpmin_;
      const FloatType stpmax_;
      int iflag_;
      bool requests_f_and_g_;
      bool requests_diag_;
      SizeType iter_;
      SizeType nfun_;
      FloatType stp_;
      FloatType stp1;
      FloatType ftol;
      FloatType ys;
      SizeType point;
      SizeType npt;
      const SizeType ispt;
      const SizeType iypt;
      int info;
      SizeType bound;
      SizeType nfev;
      std::vector<FloatType> w_;
      std::vector<FloatType> scratch_array_;
  };

  template <typename FloatType, typename SizeType>
  bool minimizer<FloatType, SizeType>::generic_run(
    FloatType* x,
    FloatType f,
    const FloatType* g,
    bool diagco,
    //---> Insertion starts
    bool gradient_only,
    bool line_search,
    //<--- Insertion ends
    const FloatType* diag)
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
        for (SizeType i = 0; i < n_; i++) {
          if (diag[i] <= FloatType(0)) {
            throw_diagonal_element_not_positive(i);
          }
        }
      }
      else {
        std::fill_n(scratch_array_.begin(), n_, FloatType(1));
        diag = &(*(scratch_array_.begin()));
      }
      for (SizeType i = 0; i < n_; i++) {
        w[ispt + i] = -g[i] * diag[i];
      }
      FloatType gnorm = std::sqrt(detail::ddot(n_, g, g));
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
        ys = detail::ddot(
          n_, w, iypt + npt, SizeType(1), w, ispt + npt, SizeType(1));
        if (!diagco) {
          FloatType yy = detail::ddot(
            n_, w, iypt + npt, SizeType(1), w, iypt + npt, SizeType(1));
          std::fill_n(scratch_array_.begin(), n_, ys / yy);
          diag = &(*(scratch_array_.begin()));
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
        if (diag == 0) {
          throw error_internal_error(__FILE__, __LINE__);
        }
        if (diagco) {
          for (SizeType i = 0; i < n_; i++) {
            if (diag[i] <= FloatType(0)) {
              throw_diagonal_element_not_positive(i);
            }
          }
        }
        SizeType cp = point;
        if (point == 0) cp = m_;
        w[n_ + cp -1] = 1 / ys;
        SizeType i;
        for (i = 0; i < n_; i++) {
          w[i] = -g[i];
        }
        cp = point;
        for (i = 0; i < bound; i++) {
          if (cp == 0) cp = m_;
          cp--;
          FloatType sq = detail::ddot(
            n_, w, ispt + cp * n_, SizeType(1), w, SizeType(0), SizeType(1));
          SizeType inmc=n_+m_+cp;
          SizeType iycn=iypt+cp*n_;
          w[inmc] = w[n_ + cp] * sq;
          detail::daxpy(n_, -w[inmc], w, iycn, w);
        }
        for (i = 0; i < n_; i++) {
          w[i] *= diag[i];
        }
        for (i = 0; i < bound; i++) {
          FloatType yr = detail::ddot(
            n_, w, iypt + cp * n_, SizeType(1), w, SizeType(0), SizeType(1));
          FloatType beta = w[n_ + cp] * yr;
          SizeType inmc=n_+m_+cp;
          beta = w[inmc] - beta;
          SizeType iscn=ispt+cp*n_;
          detail::daxpy(n_, beta, w, iscn, w);
          cp++;
          if (cp == m_) cp = 0;
        }
        std::copy(w, w+n_, w+(ispt + point * n_));
      }
      stp_ = FloatType(1);
      if (iter_ == 1) stp_ = stp1;
      std::copy(g, g+n_, w);
    }
    //---> Insertion starts
    mcsrch_instance.run(
      gradient_only, line_search, gtol_, stpmin_, stpmax_, n_, x, f, g, w, ispt + point * n_,
      stp_, ftol, xtol_, maxfev_, info, nfev, &(*(scratch_array_.begin())));
    //<--- Insertion ends
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
    for (SizeType i = 0; i < n_; i++) {
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
  template <typename FloatType, typename SizeType = std::size_t>
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
        SizeType n,
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
      SizeType n() const { return n_; }

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
        FloatType xnorm = std::sqrt(detail::ddot(n_, x, x));
        FloatType gnorm = std::sqrt(detail::ddot(n_, g, g));
        if (gnorm <= eps_ * std::max(FloatType(1), xnorm)){
          return true;
        }
        return false;
      }
    protected:
      const SizeType n_;
      const FloatType eps_;
  };

}} // namespace scitbx::lbfgs

#endif // SCITBX_LBFGS_H
