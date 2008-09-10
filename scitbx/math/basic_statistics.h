#ifndef SCITBX_MATH_BASIC_STATISTICS_H
#define SCITBX_MATH_BASIC_STATISTICS_H

#include <scitbx/array_family/misc_functions.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/math/accumulators.h>
#include <cmath>
#include <cstddef>

namespace scitbx { namespace math {

  //! Collection of basic statistics such as min, max, mean, variance.
  /*! See also:
        http://mathworld.wolfram.com/Variance.html
        http://mathworld.wolfram.com/Kurtosis.html
   */
  template <typename FloatType = double>
  class basic_statistics
  {
    public:
      typedef FloatType float_type;
      typedef std::size_t size_type;

      //! Default constructor. Some data members are not initialized!
      basic_statistics() {}

      //! Computation of all statistics.
      /*! The calculations involve two loops over the array of values.
          The memory requirements are minimal since large temporaries
          are avoided.
       */
      basic_statistics(af::const_ref<FloatType> values)
      {
        using namespace accumulator;
        n = values.size();
        if (n == 0) {
          min = -1;
          max = -1;
          max_absolute = -1;
          sum = -1;
          mean = -1;
          mean_absolute_deviation_from_mean = -1;
          biased_variance = -1;
          biased_standard_deviation = -1;
          bias_corrected_variance = -1;
          bias_corrected_standard_deviation = -1;
          skew = -1;
          kurtosis = -1;
          kurtosis_excess = -1;
          return;
        }
        min_max_accumulator<FloatType,
          mean_variance_accumulator<FloatType,
            enumerated_accumulator<FloatType> > >
              acc_1st(values[0]);
        for(std::size_t i=1; i<n; i++) acc_1st(values[i]);
        min = acc_1st.min();
        max = acc_1st.max();
        max_absolute = acc_1st.max_absolute();
        sum = acc_1st.sum();
        mean = acc_1st.mean();
        biased_variance = acc_1st.biased_variance();
        biased_standard_deviation = acc_1st.biased_standard_deviation();
        if (n == 1) {
          mean_absolute_deviation_from_mean = 0;
          bias_corrected_variance = -1;
          bias_corrected_standard_deviation = -1;
          skew = -1;
          kurtosis = -1;
          kurtosis_excess = -1;
          return;
        }
        bias_corrected_variance = acc_1st.unbiased_variance();
        bias_corrected_standard_deviation =
          acc_1st.unbiased_standard_deviation();
        if (bias_corrected_variance == 0) {
          mean_absolute_deviation_accumulator<FloatType,
            deviation_accumulator<FloatType> >
              acc_2nd(mean);
          mean_absolute_deviation_from_mean = acc_2nd.mean_absolute_deviation();
          skew = -1;
          kurtosis = -1;
          kurtosis_excess = -1;
          return;
        }
        kurtosis_accumulator<FloatType,
          skew_accumulator<FloatType,
            mean_absolute_deviation_accumulator<FloatType,
              normalised_deviation_accumulator<FloatType> > > >
                acc_2nd(mean, biased_standard_deviation);
        for(std::size_t i=0;i<n;i++) acc_2nd(values[i]);
        mean_absolute_deviation_from_mean = acc_2nd.mean_absolute_deviation();
        skew = acc_2nd.skew();
        kurtosis = acc_2nd.kurtosis();
        kurtosis_excess = acc_2nd.kurtosis_excess();
      }

      //! Number of values.
      std::size_t n;
      //! Minimum of values.
      FloatType min;
      //! Maximum of values.
      FloatType max;
      //! Maximum of absolute values.
      FloatType max_absolute;
      //! sum(values)
      FloatType sum;
      //! sum(values) / n
      FloatType mean;
      //! sum(abs(value-mean)) / n
      FloatType mean_absolute_deviation_from_mean;
      //! sum((value-mean)**2) / n
      FloatType biased_variance;
      //! sqrt(sum((value-mean)**2) / n)
      FloatType biased_standard_deviation;
      //! sum((value-mean)**2) / (n-1)
      FloatType bias_corrected_variance;
      //! sqrt(sum((value-mean)**2) / (n-1))
      FloatType bias_corrected_standard_deviation;
      //! (sum((value-mean)**3)/n) / (sum((value-mean)**2)/n)**(3/2)
      FloatType skew;
      //! (sum((value-mean)**4)/n) / (sum((value-mean)**2)/n)**2
      FloatType kurtosis;
      //! (sum((value-mean)**4)/n) / (sum((value-mean)**2)/n)**2 - 3
      FloatType kurtosis_excess;
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_BASIC_STATISTICS_H
