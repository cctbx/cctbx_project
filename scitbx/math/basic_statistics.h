#ifndef SCITBX_MATH_BASIC_STATISTICS_H
#define SCITBX_MATH_BASIC_STATISTICS_H

#include <scitbx/array_family/misc_functions.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/math/accumulators.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
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
          skewness_accumulator<FloatType,
            mean_absolute_deviation_accumulator<FloatType,
              normalised_deviation_accumulator<FloatType> > > >
                acc_2nd(mean, biased_standard_deviation);
        for(std::size_t i=0;i<n;i++) acc_2nd(values[i]);
        mean_absolute_deviation_from_mean = acc_2nd.mean_absolute_deviation();
        skew = acc_2nd.skewness();
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

  /// Median and median absolute deviation
  /** The algorithm is Quickselect, with a random choice of a pivot
      and a tweak for the case when the number of data is even,
      so as to find the both of the central values during the same search,
      the median being then the average of that pair.
  */
  template <typename FloatType = double>
  class median_statistics
  {
  public:
    typedef FloatType float_type;

    /// Construct statistics, overwriting the referenced array.
    /** It is overwritten with the absolute deviations from the median,
        but those are not stored in the same order as the original array.
    */
    median_statistics(af::ref<float_type> const &data)
      : data(data)
    {
      SCITBX_ASSERT(data.size());
      std::size_t n = data.size();

      if (n == 1) {
        median = data[0];
        median_absolute_deviation = 0;
        return;
      }

      k = data.begin() + n / 2;
      odd = n % 2;
      km1 = k-1;

      median = get_median();

      for (int i=0; i < n; ++i) {
        data[i] = std::abs(data[i] - median);
      }
      median_absolute_deviation = get_median();
    }

    /// Median of the data passed to the constructor
    float_type median;

    /// Median absolute deviation of the data passed to the constructor
    float_type median_absolute_deviation;

  private:
    af::ref<float_type> const &data;
    boost::mt19937 rnd;
    float_type *k, *km1;
    bool odd;

    float_type get_median() {
      float_type *l=data.begin(), *r=data.end()-1;

      float_type *p; // pivot
      if (odd) {
        for(;;) {
          boost::uniform_int<std::size_t> gen(0, r - l);
          p = partition(l, r, l + gen(rnd));
          if      (k < p) r = p - 1;
          else if (p < k) l = p + 1;
          else return *k;
        }
      }
      else {
        bool already_found_one_central_value = false;
        float_type central_value, other_central_value;
        for(;;) {
          boost::uniform_int<std::size_t> gen(0, r - l);
          p = partition(l, r, l + gen(rnd));
          if      (k < p)   r = p - 1;
          else if (p < km1) l = p + 1;
          else if (already_found_one_central_value) {
            other_central_value = *p;
            break;
          }
          else {
            central_value = *p;
            if   (p == k)        r = p - 1;
            else /* p == k-1 */  l = p + 1;
            already_found_one_central_value = true;
          }
        }
        return (central_value + other_central_value)/2;
      }
    }

    float_type *partition(float_type *l, float_type *r, float_type *p) {
      float_type pv = *p;
      std::swap(*p, *r);
      float_type *s = l;
      for (float_type *q=l; q < r; ++q) {
        if (*q < pv) {
          std::swap(*s, *q);
          ++s;
        }
      }
      std::swap(*r, *s);
      return s;
    }
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_BASIC_STATISTICS_H
