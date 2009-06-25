#ifndef SCITBX_LBFGSB_H
#define SCITBX_LBFGSB_H

#if defined(SCITBX_LBFGSB_HAVE_FORTRAN_LIB)
# include <scitbx/lbfgsb/fortran_interface.h>
#endif

#include <scitbx/lbfgsb/raw.h>
#include <scitbx/array_family/shared.h>

namespace scitbx {

//! C++ port of L-BFGS-B Version 2.1
/*! Original FORTRAN distribution:

      http://www.ece.northwestern.edu/~nocedal/lbfgsb.html

    Written by Ciyou Zhu in collaboration with R.H. Byrd, P. Lu-Chen
    and J. Nocedal.

    C++ port by Ralf W. Grosse-Kunstleve.
 */
namespace lbfgsb {

  //! High-level interface to L-BFGS-B Version 2.1.
  /*! An example driver can be found in the file
      scitbx/lbfgsb/boost_python/tst_lbfgsb.py,
      function driver1().
   */
  template <typename FloatType=double>
  class minimizer
  {
    public:
      //! Default constructor. Some members are not initialized!
      minimizer() {}

      //! Initialization of invariants and internal work arrays.
      /*! n is the dimension of the problem.

          @param m is the maximum number of variable metric corrections
              used to define the limited memory matrix.

          @param l is the lower bound on x (l.size() == n).

          @param u is the upper bound on x (u.size() == n).

          @param nbd represents the type of bounds imposed on the
              variables, and must be specified as follows:
                - nbd[i]=0 if x[i] is unbounded,
                - nbd[i]=1 if x[i] has only a lower bound,
                - nbd[i]=2 if x[i] has both lower and upper bounds, and
                - nbd[i]=3 if x[i] has only an upper bound.

          @param factr >= 0 is used when detecting convergence.
              The iteration will stop when

                  (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1}
                    <= factr()*floating_point_epsilon().

              Typical values for factr():
                - low accuracy:            1.e+12
                - moderate accuracy:       1.e+7
                - extremely high accuracy: 1.e+1

          @param pgtol >= 0 is used when detecting convergence.
              The iteration will stop when

                  max{|proj g[i]| i = 0, ..., n-1} <= pgtol

              where proj g[i] is the ith component of the
              projected gradient.

          @param iprint controls the frequency and type of output generated:
            iprint<0    no output is generated;
            iprint=0    print only one line at the last iteration;
            0<iprint<99 print also f and |proj g| every iprint iterations;
            iprint=99   print details of every iteration except n-vectors;
            iprint=100  print also the changes of active set and final x;
            iprint>100  print details of every iteration including x and g;
          Currently all output goes to stdout.
          When iprint > 0, the content of iterate.dat (see raw::setulb)
          is interleaved (this should be changed).

          task() is initialized as "START".
       */
      minimizer(
        int const& n,
        int const& m,
        af::shared<FloatType> l,
        af::shared<FloatType> u,
        af::shared<int> nbd,
        bool enable_stp_init,
        FloatType const& factr,
        FloatType const& pgtol,
        int const& iprint)
      :
        n_(n),
        m_(m),
        l_(l),
        u_(u),
        nbd_(nbd),
        enable_stp_init_(enable_stp_init),
        f_(0),
        factr_(factr),
        pgtol_(pgtol),
        wa_(2*m*n+4*n+12*m*m+12*m, FloatType(0)),
        iwa_(3*n, int(0)),
        task_("START"),
        iprint_(iprint),
        csave_(""),
        lsave_(4, false),
        isave_(44, int(0)),
        dsave_(29, FloatType(0)),
        requests_f_and_g_(false),
        requests_stp_init_(false),
        is_terminated_(false)
      {
        SCITBX_ASSERT(l.size() == n);
        SCITBX_ASSERT(u.size() == n);
        SCITBX_ASSERT(nbd.size() == n);
      }

      //! This routine is called repeatedly under the control of task().
      /*! Calls raw::setulb(). The return value is equivalent to
          requests_f_and_g_().
       */
      bool
      process(
        af::ref<FloatType> const& x,
        FloatType const& f,
        af::ref<FloatType> const& g,
        bool use_fortran_library=false)
      {
        SCITBX_ASSERT(!is_terminated_);
        if (task_[0] == 'F') {
          f_ = f;
          if (f_list_.size() == 0) {
            f_list_.push_back(f_);
          }
        }
        if (!use_fortran_library) {
          raw::setulb(
            n_,
            m_,
            raw::ref1<FloatType>(x),
            raw::ref1<FloatType>(l_.ref()),
            raw::ref1<FloatType>(u_.ref()),
            raw::ref1<int>(nbd_.ref()),
            f_,
            raw::ref1<FloatType>(g),
            factr_,
            pgtol_,
            raw::ref1<FloatType>(wa_.ref()),
            raw::ref1<int>(iwa_.ref()),
            task_,
            iprint_,
            csave_,
            raw::ref1<bool>(lsave_.ref()),
            raw::ref1<int>(isave_.ref()),
            raw::ref1<FloatType>(dsave_.ref()),
            enable_stp_init_);
        }
        else {
#if !defined(SCITBX_LBFGSB_HAVE_FORTRAN_LIB)
          throw error("L-BFGS-B FORTRAN library is not available.");
#else
          SCITBX_ASSERT(!enable_stp_init());
          fortran_interface::setulb(
            n_,
            m_,
            raw::ref1<FloatType>(x),
            raw::ref1<FloatType>(l_.ref()),
            raw::ref1<FloatType>(u_.ref()),
            raw::ref1<int>(nbd_.ref()),
            f_,
            raw::ref1<FloatType>(g),
            factr_,
            pgtol_,
            raw::ref1<FloatType>(wa_.ref()),
            raw::ref1<int>(iwa_.ref()),
            task_,
            iprint_,
            csave_,
            raw::ref1<bool>(lsave_.ref()),
            raw::ref1<int>(isave_.ref()),
            raw::ref1<FloatType>(dsave_.ref()));
#endif
        }
        requests_f_and_g_ = false;
        requests_stp_init_ = false;
        int t0 = task_[0];
        if (t0 == 'N') { // "NEW_X"
          f_list_.push_back(f_);
        }
        else if (t0 == 'F') { // "FG"
          requests_f_and_g_ = true;
        }
        else if (task_.substr(0,9) == "STP_INIT:") {
          requests_stp_init_ = true;
        }
        else {
          is_terminated_ = true;
          if (f_ != f_list_.back()) {
            f_list_.push_back(f_);
          }
        }
        return requests_f_and_g_;
      }

      //! Request to call process() again with newly evaluated f(x) and g(x).
      /*! Equivalent to task().substr(0,2) == "FG".
       */
      bool
      requests_f_and_g() const { return requests_f_and_g_; }

      //! Request to call process() again after optionally resetting stp.
      /*! Equivalent to task().substr(0,9) == "STP_INIT:".
       */
      bool
      requests_stp_init() const
      {
        SCITBX_ASSERT(enable_stp_init());
        return requests_stp_init_;
      }

      //! The minimization was terminated.
      /*! Use task() to obtain detailed information.
          Subsequent calls to process() will lead to an exception
          unless request_restart() was called before.
       */
      bool
      is_terminated() const { return is_terminated_; }

      //! Status of minimization process.
      /*! Possible values:

         Errors detected at start (after first call of process()):
           - ERROR: N .LE. 0
           - ERROR: M .LE. 0
           - ERROR: FACTR .LT. 0
           - ERROR: INVALID NBD
           - ERROR: NO FEASIBLE SOLUTION
           - ERROR: STP .LT. STPMIN
           - ERROR: STP .GT. STPMAX
           - ERROR: INITIAL G .GE. ZERO
           - ERROR: FTOL .LT. ZERO
           - ERROR: GTOL .LT. ZERO
           - ERROR: XTOL .LT. ZERO
           - ERROR: STPMAX .LT. STPMIN
           - ERROR: STPMIN .LT. ZERO

         User is prompted to provide f and g at current x:
           - FG_START
           - FG_LNSRCH
           - FG

         The minimizer has completed a step:
           - NEW_X

         Automatic termination (x,f,g are set to values at best solution):
           - CONVERGENCE
           - CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL
           - CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH
           - WARNING: ROUNDING ERRORS PREVENT PROGRESS
           - WARNING: STP = STPMAX
           - WARNING: STP = STPMIN
           - WARNING: XTOL TEST SATISFIED
           - ABNORMAL_TERMINATION_IN_LNSRCH

         User-defined termination via request_stop():
           - STOP: NO RESTORE

         User-defined termination via request_stop_with_restore():
           - STOP: CPU
       */
      std::string
      task() const { return task_; }

      //! List of f(x_start) and f(x) after each iteration.
      af::shared<FloatType>
      f_list() const { return f_list_; }

      //! Current value of f(x), the target function.
      /*! Last value passed to process() or as restored when terminating.
       */
      FloatType
      f() const { return f_; }

      //! Re-initialization of the minimization.
      /*! The current state is lost at the next call of process().
          f_list() is lost immediately.
       */
      void
      request_restart()
      {
        task_ = "START";
        requests_f_and_g_ = false;
        requests_stp_init_ = false;
        is_terminated_ = false;
        f_list_ = af::shared<FloatType>();
      }

      //! The minimization is stopped at the next call of process().
      /*! x(), f() and g() will not be updated.
       */
      void
      request_stop() { task_ = "STOP: NO RESTORE"; }

      //! The minimization is stopped at the next call of process().
      /*! x(), f() and g() will be restored to the values at the
          last iterate.
       */
      void
      request_stop_with_restore() { task_ = "STOP: CPU"; }

      //! Value as passed to the constructor.
      int
      n() const { return n_; }

      //! Value as passed to the constructor.
      int
      m() const { return m_; }

      //! Array of lower bounds as passed to the constructor.
      af::shared<FloatType>
      l() const { return l_; }

      //! Array of upper bounds as passed to the constructor.
      af::shared<FloatType>
      u() const { return u_; }

      //! Array of type of bounds as passed to the constructor.
      af::shared<int>
      nbd() const { return nbd_; }

      //! Value as passed to the constructor.
      bool
      enable_stp_init() const { return enable_stp_init_; }

      //! Value as passed to the constructor.
      FloatType
      factr() const { return factr_; }

      //! Value as passed to the constructor.
      FloatType
      pgtol() const { return pgtol_; }

      //! Value as passed to the constructor.
      int
      iprint() const { return iprint_; }

      /*! \brief The initial X has (not) been replaced by its projection in
          the feasible set.
       */
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      bool
      initial_x_replaced_by_projection() const { return lsave_[1-1]; }

      //! The problem is (not) constrained.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      bool
      is_constrained() const { return lsave_[2-1]; }

      //! (Not) all variables have upper and lower bounds.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      bool
      is_fully_constrained() const { return lsave_[3-1]; }

      //! The number of the current iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_iteration() { return isave_[30-1]; }

      //! Total number of function and gradient evaluations.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_fg_evaluations_total() { return isave_[34-1]; }

      //! Number of function and gradient evaluations in the current iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_fg_evaluations_iter() { return isave_[36-1]; }

      //! Total number of intervals explored in the search of Cauchy points.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_intervals_explored_cauchy_search_total() const { return isave_[22-1]; }

      /*! \brief Number of intervals explored in the search of Cauchy points
          in the current iteration.
       */
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_intervals_explored_cauchy_search_iter() const { return isave_[33-1]; }

      //! Total number of skipped BFGS updates before the current iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_skipped_bfgs_updates_total() const { return isave_[26-1]; }

      //! Total number of BFGS updates prior to the current iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_bfgs_updates_total() const { return isave_[31-1]; }

      //! Flag to indicate if subspace argmin is within or beyond the box.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      bool
      subspace_argmin_is_within_box()
      {
        if (isave_[37-1] == 0) return true;
        return false;
      }

      //! Number of free variables in the current iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_free_variables() const { return isave_[38-1]; }

      //! Number of active constraints in the current iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_active_constraints() const { return isave_[39-1]; }

      /*! \brief Number of variables leaving the set of active constraints
          in the current iteration.
       */
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_variables_leaving_active_set() const { return n_ + 1 - isave_[40-1]; }

      /*! \brief Number of variables entering the set of active constraints
          in the current iteration.
       */
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      int
      n_variables_entering_active_set() const { return isave_[41-1]; }

      //! Current 'theta' in the BFGS matrix.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      theta_bfgs_matrix_current() { return dsave_[1-1]; }

      //! f(x) in the previous iteration.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      f_previous_iteration() { return dsave_[2-1]; }

      //! Automatically determined machine precision.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
          See also: scitbx::math::floating_point_epsilon
       */
      FloatType
      floating_point_epsilon() { return dsave_[5-1]; }

      //! factr()*floating_point_epsilon().
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      factr_times_floating_point_epsilon() { return dsave_[3-1]; }

      //! 2-norm of the line search direction vector.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      two_norm_line_search_direction_vector() { return dsave_[4-1]; }

      //! Square of the 2-norm of the line search direction vector.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      two_norm_line_search_direction_vector_sq() { return dsave_[16-1]; }

      //! Accumulated time spent on searching for Cauchy points.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      accumulated_time_cauchy_search() { return dsave_[7-1]; }

      //! Accumulated time spent on subspace minimization.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      accumulated_time_subspace_minimization() { return dsave_[8-1]; }

      //! Accumulated time spent on line search.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      accumulated_time_line_search() { return dsave_[9-1]; }

      /*! \brief Slope of the line search function at the current
          point of line search.
       */
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      slope_line_search_function_current() { return dsave_[11-1]; }

      /*! \brief Slope of the line search function at the starting point
          of the line search.
       */
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      slope_line_search_function_start() { return dsave_[15-1]; }

      //! Maximum relative step length imposed in line search.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      maximum_relative_step_length() { return dsave_[12-1]; }

      //! Relative step length in the line search.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      relative_step_length_line_search() { return dsave_[14-1]; }

      //! Sets relative step length in the line search.
      /*! Should only be called if task()[:9] == "STP_INIT:"
          but this is not enforced.
       */
      void
      set_relative_step_length_line_search(FloatType const& value)
      {
        dsave_[14-1] = value;
      }

      //! Infinity norm of the projected gradient.
      /*! Should only be called if task() == "NEW_X" but this is not enforced.
       */
      FloatType
      infinity_norm_projected_gradient() { return dsave_[13-1]; }

      //! Copy of current search direction.
      af::shared<double>
      current_search_direction() // not const only b/o missing const_ref1
      {
        int ld = isave_[14-1];
        raw::ref1<FloatType> d = raw::ref1<FloatType>(wa_.ref()).get1(ld, n_);
        return af::shared<double>(d.begin(), d.end());
      }

    protected:
      int n_;
      int m_;
      af::shared<FloatType> l_;
      af::shared<FloatType> u_;
      af::shared<int> nbd_;
      bool enable_stp_init_;
      FloatType f_;
      FloatType factr_;
      FloatType pgtol_;
      af::shared<FloatType> wa_;
      af::shared<int> iwa_;
      std::string task_;
      int iprint_;
      std::string csave_;
      af::shared<bool> lsave_;
      af::shared<int> isave_;
      af::shared<FloatType> dsave_;
      bool requests_f_and_g_;
      bool requests_stp_init_;
      bool is_terminated_;
      af::shared<FloatType> f_list_;
  };

}} // namespace scitbx::lbfgsb

#endif // SCITBX_LBFGSB_H
