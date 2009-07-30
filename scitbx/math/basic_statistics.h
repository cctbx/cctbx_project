#ifndef SCITBX_MATH_BASIC_STATISTICS_H
#define SCITBX_MATH_BASIC_STATISTICS_H

#include <scitbx/array_family/misc_functions.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/math/accumulators.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <cmath>
#include <algorithm>
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


  /// Median and dispersion around it
  template <typename FloatType>
  struct median_statistics
  {
    typedef FloatType float_type;

    float_type median;

    /// The median of the absolute deviations from the median
    float_type median_absolute_deviation;

    median_statistics(float_type median, float_type median_absolute_deviation)
      : median(median), median_absolute_deviation(median_absolute_deviation)
    {}
  };


  /// functor returning a given quantile of a list of data
  /**
   For an odd number of data, the median is the central value after sorting
   the data. For an even number of data, we used the definition favoured
   by statisticians, i.e. the mean of the two central values after sorting.

   The algorithm used is Quickselect with randomisation of pivots.
   As a result, a random number generator is necessary which is part of the
   private state of instances.
   */
  class median_functor
  {
  public:
    typedef boost::mt19937 random_number_engine_t;

    /// initialise the randon number engine with a default seed
    median_functor() : rnd() {}

    /// initialise the randon number engine with the given seed
    median_functor(random_number_engine_t::result_type seed) : rnd(seed) {}

    /// the median of the given data
    template <typename FloatType>
    FloatType operator()(af::ref<FloatType> const &data) {
      SCITBX_ASSERT(data.size());
      std::size_t n = data.size();

      // early return
      if      (n == 1) return data[0];
      else if (n == 2) return (data[0] + data[1])/2;

      // hard work
      FloatType *p; // pivot

      /* for odd n, k points will point at the central value
         for even n, k will points at the farther of the two central values
       */
      FloatType * const k = data.begin() + n/2;
      FloatType *l = data.begin(), *r = data.end() - 1;
      for(;;) {
        boost::uniform_int<std::size_t> gen(0, r - l);
        p = partition(l, r, l + gen(rnd));
        if      (k < p) r = p - 1;
        else if (p < k) l = p + 1;
        else break;
      }

      // odd number of data: done
      if (n % 2) return *k;

      /* even number of data:
         the other central value is the smallest element in [data.begin(), k)
       */
      FloatType *k1 = std::max_element(data.begin(), k);
      return (*k + *k1)/2;
    }

    template <typename FloatType>
    median_statistics<FloatType>
    dispersion(af::ref<FloatType> const &data) {
      FloatType median = operator()(data);
      for (int i=0; i < data.size(); ++i) {
        data[i] = std::abs(data[i] - median);
      }
      FloatType median_absolute_deviation = operator()(data);
      return median_statistics<FloatType>(median, median_absolute_deviation);
    }

  private:
    random_number_engine_t rnd;

    /** The returned value s is obtained by re-arranging the range [l, r]
        containing p so that
        for any element x of [l, s) and for any element y of (s, r],
        *x < *s = *p < *y
     */
    template <typename FloatType>
    static
    FloatType *partition(FloatType *l, FloatType *r, FloatType *p) {
      FloatType pv = *p;
      std::swap(*p, *r);
      FloatType *s = l;
      for (FloatType *q=l; q < r; ++q) {
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
