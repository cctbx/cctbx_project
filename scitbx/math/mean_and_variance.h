#ifndef SCITBX_MATH_MEAN_AND_VARIANCE_H
#define SCITBX_MATH_MEAN_AND_VARIANCE_H

#include <scitbx/array_family/ref_reductions.h>

namespace scitbx { namespace math {

  template <typename FloatType = double>
  class mean_and_variance
  {
    public:
      typedef FloatType float_type;
      typedef std::size_t size_type;

      mean_and_variance() {}

      mean_and_variance(
        af::const_ref<FloatType> values)
      :
        have_weights_(false),
        sum_weights_(values.size()),
        sum_weights_sq_(values.size()),
        sum_weights_values_(af::sum(values)),
        sum_weights_delta_sq_(0)
      {
        FloatType m = mean();
        for(std::size_t i=0;i<values.size();i++) {
          sum_weights_delta_sq_ += fn::pow2(values[i] - m);
        }
      }

      mean_and_variance(
        af::const_ref<FloatType> values,
        af::const_ref<FloatType> weights)
      :
        have_weights_(true),
        sum_weights_(af::sum(weights)),
        sum_weights_sq_(af::sum_sq(weights)),
        sum_weights_values_(0),
        sum_weights_delta_sq_(0)
      {
        SCITBX_ASSERT(values.size() == weights.size());
        for(std::size_t i=0;i<values.size();i++) {
          sum_weights_values_ += values[i] * weights[i];
        }
        FloatType m = mean();
        for(std::size_t i=0;i<values.size();i++) {
          sum_weights_delta_sq_ += fn::pow2(values[i] - m) * weights[i];
        }
      }

      bool
      have_weights() const { return have_weights_; }

      FloatType
      mean() const
      {
        SCITBX_ASSERT(sum_weights_ > 0);
        return sum_weights_values_ / sum_weights_ ;
      }

      //! Emulation of gsl_stats_wvariance of the GNU Scientific Library.
      /*! http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
       */
      FloatType
      gsl_stats_wvariance() const
      {
        SCITBX_ASSERT(fn::pow2(sum_weights_) > sum_weights_sq_);
        return sum_weights_ / (fn::pow2(sum_weights_) - sum_weights_sq_)
             * sum_weights_delta_sq_;
      }

      //! Emulation of gsl_stats_wsd of the GNU Scientific Library.
      /*! http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html
       */
      FloatType
      gsl_stats_wsd() const { return std::sqrt(gsl_stats_wvariance()); }

      FloatType
      standard_error_of_mean_calculated_from_sample_weights() const
      {
        SCITBX_ASSERT(sum_weights_ > 0);
        return 1 / std::sqrt((sum_weights_));
      }

      FloatType
      unweighted_sample_variance() const
      {
        SCITBX_ASSERT(!have_weights_);
        SCITBX_ASSERT(sum_weights_ > 1);
        return sum_weights_delta_sq_ / (sum_weights_ - 1);
      }

      FloatType
      unweighted_sample_standard_deviation() const
      {
        return std::sqrt(unweighted_sample_variance());
      }

      FloatType
      unweighted_standard_error_of_mean() const
      {
        return std::sqrt(unweighted_sample_variance() / sum_weights_);
      }

      FloatType
      sum_weights() const { return sum_weights_; }

      FloatType
      sum_weights_sq() const { return sum_weights_sq_; }

      FloatType
      sum_weights_values() const { return sum_weights_values_; }

      FloatType
      sum_weights_delta_sq() const { return sum_weights_delta_sq_; }

    protected:
      bool have_weights_;
      FloatType sum_weights_;
      FloatType sum_weights_sq_;
      FloatType sum_weights_values_;
      FloatType sum_weights_delta_sq_;
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_MEAN_AND_VARIANCE_H
