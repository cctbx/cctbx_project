#ifndef CCTBX_MAPTBX_STATISTICS_H
#define CCTBX_MAPTBX_STATISTICS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
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

  //! Maximum Entropy Map modification (MEM): iteration update.
  inline void compute_mem_iteration (
    af::ref<double, af::c_grid_padded<3> > rho,
    af::ref<double, af::c_grid_padded<3> > delta,
    double lam,
    int n,
    double a_gd)
  {
    CCTBX_ASSERT(n>0);
    af::tiny<int, 3> n_real(af::adapt(rho.accessor().focus()));
    for (int u = 0; u < n_real[0]; u++) {
      for (int v = 0; v < n_real[1]; v++) {
        for (int w = 0; w < n_real[2]; w++) {
          double exp_arg = lam*delta(u,v,w)/n;
          double front_mul = 1.+lam/n*rho(u,v,w);
          if(delta(u,v,w)>=0) {
            double exp_ = a_gd*std::exp(-exp_arg);
            rho(u,v,w)=front_mul*exp_/(1.+lam/n*exp_);
          }
          else {
            rho(u,v,w)=front_mul*a_gd/(lam/n*a_gd+std::exp(exp_arg));
          }
        }
      }
    }
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_STATISTICS_H
