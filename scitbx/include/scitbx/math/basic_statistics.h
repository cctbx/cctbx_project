#ifndef SCITBX_MATH_BASIC_STATISTICS_H
#define SCITBX_MATH_BASIC_STATISTICS_H

#include <scitbx/array_family/misc_functions.h>
#include <cmath>

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
        }
        else {
          min = max = sum = values[0];
          for(std::size_t i=1;i<n;i++) {
            FloatType const& v = values[i];
            if (min > v) min = v;
            if (max < v) max = v;
            sum += v;
          }
          if (-min > max) max_absolute = -min;
          else            max_absolute = max;
          mean = sum / static_cast<FloatType>(n);
          FloatType sum_a = 0;
          FloatType sum_2 = 0;
          FloatType sum_3 = 0;
          FloatType sum_4 = 0;
          for(std::size_t i=0;i<n;i++) {
            FloatType const& v = values[i];
            FloatType vm = v - mean;
            FloatType vms = vm * vm;
            sum_a += fn::absolute(vm);
            sum_2 += vms;
            sum_3 += vms * vm;
            sum_4 += vms * vms;
          }
          FloatType nf = static_cast<FloatType>(n);
          mean_absolute_deviation_from_mean = sum_a / nf;
          biased_variance = sum_2 / nf;
          biased_standard_deviation = std::sqrt(biased_variance);
          if (n == 1) {
            bias_corrected_variance = -1;
            bias_corrected_standard_deviation = -1;
            skew = -1;
            kurtosis = -1;
            kurtosis_excess = -1;
          }
          else {
            bias_corrected_variance = sum_2 / static_cast<FloatType>(n-1);
            bias_corrected_standard_deviation = std::sqrt(
              bias_corrected_variance);
            skew = (sum_3 / nf) / fn::pow3(biased_standard_deviation);
            kurtosis = (sum_4 / nf) / fn::pow2(biased_variance);
            kurtosis_excess = kurtosis - 3;
          }
        }
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
