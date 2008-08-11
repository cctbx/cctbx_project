#ifndef SCITBX_LBFGS_DROP_CONVERGENCE_TEST_H
#define SCITBX_LBFGS_DROP_CONVERGENCE_TEST_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/math/linear_regression.h>

namespace scitbx { namespace lbfgs {

  //! Convergence test based on monitoring an objective function.
  /*! This test monitors an objective function <code>f</code>
      in the course of the minimization to determine when it
      has reached a plateau. The test is independent of the
      absolute scale of the objective function.

      <ul>
      <li>The maximum value of the drop in the objective function
          is determined (max_drop()).
          <p>
      <li>A straight line is fitted to the last n_test_points() values
          of the objective function. Let <code>slope</code>
          be the slope of this fitted line. Convergence is
          detected if:
          <pre>
          -slope &lt;=   max_drop() * max_drop_eps()
                       * number_of_iterations^iteration_coefficient()
          </pre>
          Note that with iteration_coefficient() > 1 this
          test becomes increasingly tolerant with the
          number of iterations. The rationale is that the
          slope will not change very much for parameters
          that are already close to the optimum.
      </ul>
   */
  template <typename FloatType, typename SizeType = std::size_t>
  class drop_convergence_test
  {
    public:
      //! Constructor.
      /*! See the class details for an outline of the convergence
          detection algorithm and the meaning of the parameters.
       */
      explicit
      drop_convergence_test(
        SizeType n_test_points=5,
        FloatType max_drop_eps=FloatType(1.e-5),
        FloatType iteration_coefficient=FloatType(2))
      : n_test_points_(n_test_points),
        max_drop_eps_(max_drop_eps),
        iteration_coefficient_(iteration_coefficient),
        max_drop_(0),
        max_f_(0)
      {
        SCITBX_ASSERT(n_test_points >= 2);
        SCITBX_ASSERT(max_drop_eps_ >= FloatType(0));
        SCITBX_ASSERT(iteration_coefficient_ >= FloatType(1));
      }

      /*! \brief Number of most recent objective function values used
           in fit of straight lines (as passed to the constructor).
       */
      SizeType n_test_points() const { return n_test_points_; }

      /*! \brief Base tolerance for test of drop of objective function
           (as passed to the constructor).
       */
      FloatType max_drop_eps() const { return max_drop_eps_; }

      /*! \brief Tolerance for comparison of slopes
           (as passed to the constructor).
       */
      FloatType iteration_coefficient() const {
        return iteration_coefficient_;
      }

      //! Execution of the convergence test.
      /*! Note that at least n_test_points() executions are required before
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
      const SizeType n_test_points_;
      const FloatType max_drop_eps_;
      const FloatType iteration_coefficient_;
      af::shared<FloatType> x_;
      af::shared<FloatType> y_;
      FloatType max_drop_;
      FloatType max_f_;
  };

  template <typename FloatType, typename SizeType>
  bool
  drop_convergence_test<FloatType, SizeType>::operator()(FloatType f)
  {
    if (x_.size()) {
      max_drop_ = std::max(max_drop_, y_.back() - f);
    }
    max_f_ = std::max(max_f_, fn::absolute(f));
    x_.push_back(x_.size() + 1);
    y_.push_back(f);
    if (x_.size() < n_test_points_) return false;
    if (!max_f_) return true; // y_ must be all 0
    af::shared<FloatType> y_scaled;
    y_scaled.reserve(n_test_points_);
    for(std::size_t i=y_.size()-n_test_points_;i<y_.size();i++) {
      y_scaled.push_back(y_[i] / max_f_);
    }
    // fit last n_test_points_ to straight line: y = m*x + b
    math::linear_regression<FloatType> linreg_y(
      af::const_ref<FloatType>(x_.end() - n_test_points_, n_test_points_),
      af::const_ref<FloatType>(y_scaled.begin(), n_test_points_));
    SCITBX_ASSERT(linreg_y.is_well_defined());
    // check absolute value of slope
    FloatType sliding_tolerance =
      max_drop_ * max_drop_eps_ * std::pow(
        FloatType(x_.size()), iteration_coefficient_);
    if (-linreg_y.slope() * max_f_ <= sliding_tolerance) {
      return true;
    }
    return false;
  }

}} // namespace scitbx::lbfgs

#endif // SCITBX_LBFGS_DROP_CONVERGENCE_TEST_H
