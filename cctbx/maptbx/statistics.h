#ifndef CCTBX_MAPTBX_STATISTICS_H
#define CCTBX_MAPTBX_STATISTICS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/accumulators.h>
#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace maptbx {

  namespace details {
    template <class AccumulatorType>
    struct generic_statistics
    {
      generic_statistics() : accumulator(0) {}

      template <typename FloatType>
      generic_statistics(af::const_ref<FloatType, af::flex_grid<> > const& map)
        : accumulator(0)
      {
        CCTBX_ASSERT(map.accessor().focus_size_1d() > 0);
        if (!map.accessor().is_padded()) {
          accumulator = AccumulatorType(map[0]);
          for(std::size_t i=1; i < map.size(); ++i) accumulator(map[i]);
        }
        else {
          typedef typename af::flex_grid<>::index_type index_type;
          af::flex_grid<> zero_based = map.accessor().shift_origin();
          af::nested_loop<index_type> iter(zero_based.focus());
          accumulator = AccumulatorType(map[zero_based(iter())]);
          while (iter.incr()) accumulator(map[zero_based(iter())]);
        }
      }

      AccumulatorType accumulator;
    };

    template <class AccumulatorType>
    struct generic_statistical_moments
    {
      generic_statistical_moments() : accumulator(0,1) {}

      template <typename FloatType, typename OtherFloatType>
      generic_statistical_moments(
        af::const_ref<OtherFloatType, af::flex_grid<> > const& map,
        FloatType about, FloatType width
      )
        : accumulator(about, width)
      {
        CCTBX_ASSERT(map.accessor().focus_size_1d() > 0);
        if (width == 0) return;
        if (!map.accessor().is_padded()) {
          for(std::size_t i=0; i < map.size(); ++i) accumulator(map[i]);
        }
        else {
          typedef typename af::flex_grid<>::index_type index_type;
          af::flex_grid<> zero_based = map.accessor().shift_origin();
          for(af::nested_loop<index_type> iter(zero_based.focus());
              !iter.over(); iter.incr()) accumulator(map[zero_based(iter())]);
        }
      }

      AccumulatorType accumulator;
    };

    using namespace scitbx::math::accumulator;

    template <typename FloatType = double>
    struct statistics_traits {
      typedef generic_statistics<
                min_max_accumulator<FloatType,
                  mean_variance_accumulator<FloatType,
                    enumerated_accumulator<FloatType> > > >
              basic_statistics_t;
      typedef generic_statistical_moments<
                skewness_accumulator<FloatType,
                  kurtosis_accumulator<FloatType,
                    normalised_deviation_accumulator<FloatType> > > >
              extra_statistics_t;
    };
  }

  //! Determines simple map statistics.
  template <typename FloatType = double>
  class statistics
  {
    public:
      //! Default constructor. Data members are not initialized!
      statistics() : basic_stats() {}

      //! Computes the statistics.
      template <typename OtherFloatType>
      statistics(af::const_ref<OtherFloatType, af::flex_grid<> > const& map)
        : basic_stats(map)
      {}

      FloatType
      min() const { return basic_stats.accumulator.min(); }

      FloatType
      max() const { return basic_stats.accumulator.max(); }

      FloatType
      mean() const { return basic_stats.accumulator.mean(); }

      FloatType
      mean_sq() const { return basic_stats.accumulator.mean_squares(); }

      FloatType
      sigma() const {
        return basic_stats.accumulator.biased_standard_deviation();
      }

    private:
      typedef typename details::statistics_traits<FloatType>::basic_statistics_t
              basic_statistics_t;
      basic_statistics_t basic_stats;
  };


  //! Determines higher order statistical central moments.
  template <typename FloatType = double>
  class more_statistics : public statistics<FloatType>
  {
    public:
      //! Default constructor. Data members are not initialized!
      more_statistics() : statistics<FloatType>(), extra_stats() {}

      //! Computes the statistics.
      template <typename OtherFloatType>
      more_statistics(af::const_ref<OtherFloatType, af::flex_grid<> > const& map)
        : statistics<FloatType>(map),
          extra_stats(map, this->mean(), this->sigma())
      {}

      FloatType
      skewness() const { return extra_stats.accumulator.skewness(); }

      FloatType
      kurtosis() const { return extra_stats.accumulator.kurtosis(); }

    private:
      typedef typename details::statistics_traits<FloatType>::extra_statistics_t
              extra_statistics_t;
      extra_statistics_t extra_stats;
  };

  // maxent stuff
  inline void normalize_and_combine (
    af::versa<double,  af::flex_grid<> > priorA_map,
    af::const_ref<double,  af::flex_grid<> > priorB_map,
    double norm,
    double current_lambda)
  {
    for (int i = 0; i < priorA_map.size(); i++) {
      double val = priorA_map[i] * norm;
      priorA_map[i] = priorB_map[i] * exp(current_lambda * val);
    }
  }

  inline double calculate_entropy (
    af::const_ref<double,  af::flex_grid<> > const& map_data)
  {
    af::tiny<int, 3> n_real(af::adapt(map_data.accessor().focus()));
    double sum = 0.0;
    double ent = 0.0;
    for (int u = 0; u < n_real[0]; u++) {
      for (int v = 0; v < n_real[1]; v++) {
        for (int w = 0; w < n_real[2]; w++) {
          sum += map_data(u,v,w);
    }}}
    for (int u = 0; u < n_real[0]; u++) {
      for (int v = 0; v < n_real[1]; v++) {
        for (int w = 0; w < n_real[2]; w++) {
          double val = map_data(u,v,w);
          ent += (-val / sum) * log(val / sum);
    }}}
    return ent;
  }

  class update_prior
  {
    public :
      double sum;
      double chi2;

      update_prior (
        af::const_ref<std::complex<double>, af::flex_grid<> > const& fobs,
        af::const_ref<std::complex<double>, af::flex_grid<> > const& sigf,
        af::versa<std::complex<double>, af::flex_grid<> > priorA)
      {
        sum = 0.0;
        chi2 = 0.0;
        for (int i = 0; i < sigf.size(); i++) {
          double sigma = sigf[i].real();
          if (sigma != 0.0) {
            double s = std::pow(sigma, 2);
            double ac = priorA[i].real();
            double bc = priorA[i].imag();
            double ao = fobs[i].real();
            double bo = fobs[i].imag();
            double amp1 = sqrt(ac*ac + bc*bc);
            double amp2 = sqrt(ao*ao + bo*bo);
            chi2 += ((amp1 - amp2) * (amp1 - amp2)) / s;
            sum += abs(abs(amp2) - abs(amp1));
            priorA[i] = std::complex<double>((ao-ac)/s, (bo-bc)/s);
          } else {
            priorA[i] = std::complex<double>(0.0,0.0);
          }
        }
        return;
      }
  };

  inline void clear_map (
    af::versa<double, af::flex_grid<> > map_data,
    double mean_density)
  {
    for (int i = 0; i < map_data.size(); i++) {
      map_data[i] = mean_density;
    }
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_STATISTICS_H
