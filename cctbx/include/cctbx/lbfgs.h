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
    Ralf W. Grosse-Kunstleve: C++ port, March 2002.
 */

#include <exception>
#include <vector>
#include <string>
#include <stdio.h>
#include <boost/config.hpp>
#include <cctbx/fixes/cstdlib>
#include <cctbx/fixes/cmath>

namespace cctbx {

  class lbfgs_error : public std::exception {
    public:
      lbfgs_error(const std::string& msg) throw()
        : msg_("lbfgs error: " + msg)
      {}
      virtual ~lbfgs_error() throw() {}
      virtual const char* what() const throw() { return msg_.c_str(); }
    protected:
      std::string msg_;
    public:
      static std::string itoa(unsigned long i) {
        char buf[80];
        sprintf(buf, "%lu", i); // FUTURE: use C++ facility
        return std::string(buf);
      }
  };

  class lbfgs_error_internal_error : public lbfgs_error {
    public:
      lbfgs_error_internal_error(const char* file, unsigned long line) throw()
        : lbfgs_error(
            "Internal Error: " + std::string(file) + "(" + itoa(line) + ")")
      {}
  };

  class lbfgs_error_improper_input_parameter : public lbfgs_error {
    public:
      lbfgs_error_improper_input_parameter(const std::string& msg) throw()
        : lbfgs_error("Improper input parameter: " + msg)
      {}
  };

  class lbfgs_error_improper_input_data : public lbfgs_error {
    public:
      lbfgs_error_improper_input_data(const std::string& msg) throw()
        : lbfgs_error("Improper input data: " + msg)
      {}
  };

  class lbfgs_error_search_direction_not_descent : public lbfgs_error {
    public:
      lbfgs_error_search_direction_not_descent() throw()
        : lbfgs_error("The search direction is not a descent direction.")
      {}
  };

  class lbfgs_error_line_search_failed : public lbfgs_error {
    public:
      lbfgs_error_line_search_failed(const std::string& msg) throw()
        : lbfgs_error("Line search failed: " + msg)
      {}
  };

  /*! This class solves the unconstrained minimization problem
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
      <p>
      The user is required to calculate the function value
      <code>f</code> and its gradient <code>g</code>. In order to
      allow the user complete control over these computations,
      reverse communication is used. The routine must be called
      repeatedly under the control of the parameter
      <code>iflag</code>.
      <p>
      The steplength is determined at each iteration by means of
      the line search routine <code>mcsrch</code>, which is a
      slight modification of the routine <code>CSRCH</code> written
      by More' and Thuente.
      <p>
      The only variables that are machine-dependent are
      <code>xtol</code>, <code>stpmin</code> and
      <code>stpmax</code>.
      <p>
      Fatal errors cause exceptions to be thrown.
   */
  class lbfgs
  {
    public:
      //! Default constructor. Members are not initialized!
      lbfgs()
      : n_(0), m_(0), eps_(0), xtol_(0), maxfev_(0),
        gtol_(0), stpmin_(0), stpmax_(0),
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

          @param maxfev Termination occurs when the number of evaluations
             of the objective function is at least <code>maxfev</code> by
             the end of an iteration.

          <code>gtol</code> controls the accuracy of the line search
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
        std::size_t n,
        std::size_t m = 5,
        double eps = 1.e-5,
        double xtol = 1.e-16,
        int maxfev = 20,
        double gtol = 0.9,
        double stpmin = 1.e-20,
        double stpmax = 1.e20)
        : n_(n), m_(m), eps_(eps), xtol_(xtol), maxfev_(maxfev),
          gtol_(gtol), stpmin_(stpmin), stpmax_(stpmax),
          iflag_(0), iter_(0), nfun_(0), gnorm_(0), stp_(0),
          stp1(0), ftol(0.0001), ys(0), point(0), npt(0),
          ispt(n+2*m), iypt((n+2*m)+n*m),
          info(0), bound(0), nfev(0)
      {
        if (n_ == 0) {
          throw lbfgs_error_improper_input_parameter("n is not positive.");
        }
        if (m_ == 0) {
          throw lbfgs_error_improper_input_parameter("m is not positive.");
        }
        if (gtol_ <= 1.e-4) {
          throw lbfgs_error_improper_input_parameter("gtol <= 1.e-4.");
        }
        w_.resize(n_*(2*m_+1)+2*m_);
      }

      //! Number of free parameters.
      std::size_t n() const { return n_; }

      //! Number of corrections kept.
      std::size_t m() const { return m_; }

      //! XXX
      double eps() const { return eps_; }

      //! XXX
      double xtol() const { return xtol_; }

      //! XXX
      int maxfev() const { return maxfev_; }

      //! XXX
      double gtol() const { return gtol_; }

      //! XXX
      double stpmin() const { return stpmin_; }

      //! XXX
      double stpmax() const { return stpmax_; }

      //! XXX
      bool is_converged() const { return iflag_ == 0 && nfun_ > 0; }

      //! XXX
      bool requests_f_and_g() const { return iflag_ == 1; }

      //! XXX
      bool requests_diag() const { return iflag_ == 2; }

      //! Number of iterations so far.
      std::size_t iter() const { return iter_; }

      //! Number of function evaluations so far.
      /*! This method returns the total number of evaluations of the
          objective function. The total number of function evaluations
          increases by the number of evaluations required for the line
          search; the total is only increased after a successful line
          search.
        */
      std::size_t nfun() const { return nfun_; }

      //! Norm of gradient at current solution <code>x</code>.
      double gnorm() const { return gnorm_; }

      //! Norm of gradient given gradient array of length n().
      double gnorm(const double* g) const {
        return std::sqrt(ddot(n_, g, 0, 1, g, 0, 1));
      }

      //! Current stepsize.
      double stp() const { return stp_; }

      //! Execution of one step of the minimization.
      /*! @param x On initial entry this must be set by the user to the
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

          <code>lbfgs</code> will use a default value for diag
          as described below. XXX where?
       */
      void run(
        double* x,
        double f,
        const double* g)
      {
        if (diag_.size() == 0) diag_.resize(n_);
        generic_run(x, f, g, false, &(*(diag_.begin())));
      }

      //! Execution of one step of the minimization.
      /*! @param x See other overload.

          @param f See other overload.

          @param g See other overload.

          @param diag On initial entry or on re-entry with
             <code>iflag = 2</code>, <code>diag</code> must be set by
             the user to contain the values of the diagonal matrix
             <code>Hk0</code>.
             The routine will return at each iteration of the algorithm
             with <code>iflag = 2</code>.
             Restriction: all elements of <code>diag</code> must be
             positive.
       */
      void run(
        double* x,
        double f,
        const double* g,
        double* diag)
      {
        generic_run(x, f, g, true, diag);
      }

    protected:
      static void throw_diagonal_element_not_positive(std::size_t i) {
        throw lbfgs_error_improper_input_data(
          "The " + lbfgs_error::itoa(i) + ". diagonal element of the"
          " inverse Hessian approximation is not positive.");
      }

      void generic_run(
        double* x,
        double f,
        const double* g,
        bool diagco,
        double* diag)
      {
        double* w = &(*(w_.begin()));
        bool execute_entire_while_loop = false;
        if (iflag_ == 0) { // Initialize.
          nfun_ = 1;
          if (diagco) {
            for (std::size_t i = 0; i < n_; i++) {
              if (diag[i] <= 0) {
                throw_diagonal_element_not_positive(i);
              }
            }
          }
          else {
            std::fill_n(diag, n_, double(1));
          }
          for (std::size_t i = 0; i < n_; i++) {
            w[ispt + i] = -g[i] * diag[i];
          }
          gnorm_ = std::sqrt(ddot(n_, g, 0, 1, g, 0, 1));
          stp1= 1/gnorm_;
          execute_entire_while_loop = true;
        }
        for (;;) {
          if (execute_entire_while_loop) {
            bound = iter_;
            iter_++;
            info = 0;
            if (iter_ != 1) {
              if (iter_ > m_) bound = m_;
              ys = ddot(n_, w, iypt + npt, 1, w, ispt + npt, 1);
              if (!diagco) {
                double yy = ddot(n_, w, iypt + npt, 1, w, iypt + npt, 1);
                std::fill_n(diag, n_, ys / yy);
              }
              else {
                iflag_ = 2;
                return;
              }
            }
          }
          if (execute_entire_while_loop || iflag_ == 2) {
            if (iter_ != 1) {
              if (diagco) {
                for (std::size_t i = 0; i < n_; i++) {
                  if (diag[i] <= 0) {
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
                double sq = ddot(n_, w, ispt + cp * n_, 1, w, 0, 1);
                std::size_t inmc=n_+m_+cp;
                std::size_t iycn=iypt+cp*n_;
                w[inmc] = w[n_ + cp] * sq;
                daxpy(n_, -w[inmc], w, iycn, 1, w, 0, 1);
              }
              for (i = 0; i < n_; i++) {
                w[i] *= diag[i];
              }
              for (i = 0; i < bound; i++) {
                double yr = ddot(n_, w, iypt + cp * n_, 1, w, 0, 1);
                double beta = w[n_ + cp] * yr;
                std::size_t inmc=n_+m_+cp;
                beta = w[inmc] - beta;
                std::size_t iscn=ispt+cp*n_;
                daxpy(n_, beta, w, iscn, 1, w, 0, 1);
                cp++;
                if (cp == m_) cp = 0;
              }
              std::copy(w, w+n_, w+(ispt + point * n_));
            }
            stp_ = double(1);
            if (iter_ == 1) stp_ = stp1;
            std::copy(g, g+n_, w);
          }
          mcsrch_instance.run(
            gtol_, stpmin_, stpmax_, n_, x, f, g, w, ispt + point * n_,
            stp_, ftol, xtol_, maxfev_, info, nfev, diag);
          if (info == -1) {
            iflag_ = 1;
            return;
          }
          if (info != 1) {
            throw lbfgs_error_internal_error(__FILE__, __LINE__);
          }
          nfun_ += nfev;
          npt = point*n_;
          for (std::size_t i = 0; i < n_; i++) {
            w[ispt + npt + i] = stp_ * w[ispt + npt + i];
            w[iypt + npt + i] = g[i] - w[i];
          }
          point++;
          if (point == m_) point = 0;
          gnorm_ = std::sqrt(ddot(n_, g, 0, 1, g, 0, 1));
          double xnorm = std::max(
            double(1), std::sqrt(ddot(n_, x, 0, 1, x, 0, 1)));
          if (gnorm_ / xnorm <= eps_) {
            iflag_ = 0;
            break;
          }
          execute_entire_while_loop = true; // from now on, execute whole loop
        }
      }

      /* Compute the sum of a vector times a scalara plus another vector.
         Adapted from the subroutine <code>daxpy</code> in
         <code>lbfgs.f</code>.
       */
      static void daxpy(
        std::size_t n,
        double da,
        const double* dx,
        std::size_t ix0,
        std::size_t incx,
        double* dy,
        std::size_t iy0,
        std::size_t incy)
      {
        std::size_t i, ix, iy, m;
        if (n == 0) return;
        if (da == double(0)) return;
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
      static double ddot(
        std::size_t n,
        const double* dx,
        std::size_t ix0,
        std::size_t incx,
        const double* dy,
        std::size_t iy0,
        std::size_t incy)
      {
        std::size_t i, ix, iy, m, mp1;
        double dtemp(0);
        if (n == 0) return double(0);
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

      // This class implements an algorithm for multi-dimensional line search.
      class mcsrch
      {
        protected:
          int infoc;
          double dginit;
          bool brackt;
          bool stage1;
          double finit;
          double dgtest;
          double width;
          double width1;
          double stx;
          double fx;
          double dgx;
          double sty;
          double fy;
          double dgy;
          double stmin;
          double stmax;

          static double sqr(double x) { return x * x; }

          static double max3(double x, double y, double z) {
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
            const double& gtol,
            const double& stpmin,
            const double& stpmax,
            std::size_t n,
            double* x,
            double f,
            const double* g,
            double* s,
            std::size_t is0,
            double& stp,
            double ftol,
            double xtol,
            int maxfev,
            int& info,
            std::size_t& nfev,
            double* wa)
          {
            if (info != -1) {
              infoc = 1;
              if (   n <= 0
                  || xtol < 0
                  || maxfev <= 0
                  || gtol < 0
                  || stpmin < 0
                  || stpmax < stpmin) {
                throw lbfgs_error_internal_error(__FILE__, __LINE__);
              }
              if (stp <= 0 || ftol < 0) {
                throw lbfgs_error_internal_error(__FILE__, __LINE__);
              }
              // Compute the initial gradient in the search direction
              // and check that s is a descent direction.
              dginit = 0;
              for (std::size_t j = 0; j < n; j++) {
                dginit += g[j] * s[is0+j];
              }
              if (dginit >= 0) {
                throw lbfgs_error_search_direction_not_descent();
              }
              brackt = false;
              stage1 = true;
              nfev = 0;
              finit = f;
              dgtest = ftol*dginit;
              width = stpmax - stpmin;
              width1 = double(2) * width;
              std::copy(x, x+n, wa);
              // The variables stx, fx, dgx contain the values of the step,
              // function, and directional derivative at the best step.
              // The variables sty, fy, dgy contain the value of the step,
              // function, and derivative at the other endpoint of
              // the interval of uncertainty.
              // The variables stp, f, dg contain the values of the step,
              // function, and derivative at the current step.
              stx = 0;
              fx = finit;
              dgx = dginit;
              sty = 0;
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
                  stmax = stp + double(4) * (stp - stx);
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
              double dg(0);
              for (std::size_t j = 0; j < n; j++) {
                dg += g[j] * s[is0+j];
              }
              double ftest1 = finit + stp*dgtest;
              // Test for convergence.
              if ((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0) {
                throw lbfgs_error_line_search_failed(
                  "Rounding errors prevent further progress."
                  " There may not be a step which satisfies the"
                  " sufficient decrease and curvature conditions."
                  " Tolerances may be too small.");
              }
              if (stp == stpmax && f <= ftest1 && dg <= dgtest) {
                throw lbfgs_error_line_search_failed(
                  "The step is at the upper bound stpmax().");
              }
              if (stp == stpmin && (f > ftest1 || dg >= dgtest)) {
                throw lbfgs_error_line_search_failed(
                  "The step is at the lower bound stpmin().");
              }
              if (nfev >= maxfev) {
                throw lbfgs_error_line_search_failed(
                  "Number of function evaluations has reached maxfev().");
              }
              if (brackt && stmax - stmin <= xtol * stmax) {
                throw lbfgs_error_line_search_failed(
                  "Relative width of the interval of uncertainty"
                  " is at most xtol().");
              }
              // Check for termination.
              if (f <= ftest1 && std::fabs(dg) <= gtol * (-dginit)) {
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
                double fm = f - stp*dgtest;
                double fxm = fx - stx*dgtest;
                double fym = fy - sty*dgtest;
                double dgm = dg - dgtest;
                double dgxm = dgx - dgtest;
                double dgym = dgy - dgtest;
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
                if (std::fabs(sty - stx) >= 0.66 * width1) {
                  stp = stx + double(0.5) * (sty - stx);
                }
                width1 = width;
                width = std::fabs(sty - stx);
              }
            }
          }

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
            double& stx,
            double& fx,
            double& dx,
            double& sty,
            double& fy,
            double& dy,
            double& stp,
            double fp,
            double dp,
            bool& brackt,
            double stpmin,
            double stpmax)
          {
            bool bound;
            double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
            int info = 0;
            if (   (   brackt && (stp <= std::min(stx, sty)
                    || stp >= std::max(stx, sty)))
                || dx * (stp - stx) >= 0.0 || stpmax < stpmin) {
              return 0;
            }
            // Determine if the derivatives have opposite sign.
            sgnd = dp * (dx / std::fabs(dx));
            if (fp > fx) {
              // First case. A higher function value.
              // The minimum is bracketed. If the cubic step is closer
              // to stx than the quadratic step, the cubic step is taken,
              // else the average of the cubic and quadratic steps is taken.
              info = 1;
              bound = true;
              theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
              s = max3(std::fabs(theta), std::fabs(dx), std::fabs(dp));
              gamma = s * std::sqrt(sqr(theta / s) - (dx / s) * (dp / s));
              if (stp < stx) gamma = - gamma;
              p = (gamma - dx) + theta;
              q = ((gamma - dx) + gamma) + dp;
              r = p/q;
              stpc = stx + r * (stp - stx);
              stpq = stx
                + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2)
                  * (stp - stx);
              if (std::fabs(stpc - stx) < std::fabs(stpq - stx)) {
                stpf = stpc;
              }
              else {
                stpf = stpc + (stpq - stpc) / 2;
              }
              brackt = true;
            }
            else if (sgnd < 0.0) {
              // Second case. A lower function value and derivatives of
              // opposite sign. The minimum is bracketed. If the cubic
              // step is closer to stx than the quadratic (secant) step,
              // the cubic step is taken, else the quadratic step is taken.
              info = 2;
              bound = false;
              theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
              s = max3(std::fabs(theta), std::fabs(dx), std::fabs(dp));
              gamma = s * std::sqrt(sqr(theta / s) - (dx / s) * (dp / s));
              if (stp > stx) gamma = - gamma;
              p = (gamma - dp) + theta;
              q = ((gamma - dp) + gamma) + dx;
              r = p/q;
              stpc = stp + r * (stx - stp);
              stpq = stp + (dp / (dp - dx)) * (stx - stp);
              if (std::fabs(stpc - stp) > std::fabs(stpq - stp)) {
                stpf = stpc;
              }
              else {
                stpf = stpq;
              }
              brackt = true;
            }
            else if (std::fabs(dp) < std::fabs(dx)) {
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
              theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
              s = max3(std::fabs(theta), std::fabs(dx), std::fabs(dp));
              gamma = s * std::sqrt(
                std::max(double(0), sqr(theta / s) - (dx / s) * (dp / s)));
              if (stp > stx) gamma = -gamma;
              p = (gamma - dp) + theta;
              q = (gamma + (dx - dp)) + gamma;
              r = p/q;
              if (r < 0.0 && gamma != 0.0) {
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
                if (std::fabs(stp - stpc) < std::fabs(stp - stpq)) {
                  stpf = stpc;
                }
                else {
                  stpf = stpq;
                }
              }
              else {
                if (std::fabs(stp - stpc) > std::fabs(stp - stpq)) {
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
                theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
                s = max3(std::fabs(theta), std::fabs(dy), std::fabs(dp));
                gamma = s * std::sqrt(sqr(theta / s) - (dy / s) * (dp / s));
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
              if (sgnd < 0.0) {
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
                stp = std::min(stx + 0.66 * (sty - stx), stp);
              }
              else {
                stp = std::max(stx + 0.66 * (sty - stx), stp);
              }
            }
            return info;
          }
      };

      mcsrch mcsrch_instance;
      const std::size_t n_;
      const std::size_t m_;
      const double eps_;
      const double xtol_;
      const int maxfev_;
      const double gtol_;
      const double stpmin_;
      const double stpmax_;
      int iflag_;
      std::size_t iter_;
      std::size_t nfun_;
      double gnorm_;
      double stp_;
      double stp1;
      double ftol;
      double ys;
      std::size_t point;
      std::size_t npt;
      const std::size_t ispt;
      const std::size_t iypt;
      int info;
      std::size_t bound;
      std::size_t nfev;
      std::vector<double> w_;
      std::vector<double> diag_;
  };

} // namespace cctbx

#endif // CCTBX_LBFGS_H
