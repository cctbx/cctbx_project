#ifndef CCTBX_MAPTBX_STATISTICS_H
#define CCTBX_MAPTBX_STATISTICS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/c_grid.h>
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

  //! Compute MEM iteration, total charge of rho_tilda (Z), total power (TP),
  //! working (Hw) and normalized (Hn) entropy.
  template <typename FloatType = double>
  class mem_iteration
  {
    FloatType scale_, tp_, z_, hw_, hn_;

    public:
      mem_iteration(
        af::ref<FloatType, af::c_grid<3> > const& rho_mod,
        af::ref<FloatType, af::c_grid<3> > const& rho_obs,
        af::ref<FloatType, af::c_grid<3> > rho,
        FloatType lam,
        af::tiny<int, 3> const& n_real,
        FloatType a_gd,
        FloatType beta,
        bool use_scale)
      :
      scale_(1), tp_(0), z_(0), hw_(0), hn_(0)
      {
        CCTBX_ASSERT(rho_mod.size() == rho_obs.size());
        CCTBX_ASSERT(rho_mod.size() == rho.size());
        if(use_scale) {
          FloatType num=0, den=0;
          for (int u = 0; u < n_real[0]; u++) {
            for (int v = 0; v < n_real[1]; v++) {
              for (int w = 0; w < n_real[2]; w++) {
                FloatType rm = std::abs(rho_mod(u,v,w));
                FloatType ro = std::abs(rho_obs(u,v,w));
                num += (rm*ro);
                den += (ro*ro);
          }}}
          if(den!=0 && num!=0) scale_ = 1./(num/den);
        }
        int cntr_rho_positive = 0;
        FloatType rho_positive_sum = 0;
        for (int u = 0; u < n_real[0]; u++) {
          for (int v = 0; v < n_real[1]; v++) {
            for (int w = 0; w < n_real[2]; w++) {
              FloatType delta = rho_mod(u,v,w) - rho_obs(u,v,w)*scale_;
              FloatType rho_tilda;
              FloatType exp_arg = lam*delta;
              FloatType front_mul = 1.+lam*rho(u,v,w);
              if(delta>=0) {
                FloatType exp_ = a_gd*std::exp(-exp_arg);
                rho_tilda=front_mul*exp_/(1.+lam*exp_);
              }
              else rho_tilda=front_mul*a_gd/(lam*a_gd+std::exp(exp_arg));
              z_ += rho_tilda;
              FloatType rho_new = (1-beta)*rho(u,v,w) + beta*rho_tilda;
              rho(u,v,w) = rho_new;
              tp_ += rho_new;
              if(rho_new > 0) {
                hw_ += rho_new * std::log(rho_new);
                cntr_rho_positive += 1;
                rho_positive_sum += rho_new;
        }}}}
        hw_ = -hw_;
        for (int u = 0; u < n_real[0]; u++) {
          for (int v = 0; v < n_real[1]; v++) {
            for (int w = 0; w < n_real[2]; w++) {
              FloatType rho_uvw = rho(u,v,w);
              if(rho_uvw > 0) {
                FloatType rho_over_sum_rho = rho_uvw/rho_positive_sum;
                hn_ += rho_over_sum_rho * std::log(rho_over_sum_rho);
        }}}}
        hn_ = -hn_/std::log((FloatType) cntr_rho_positive);
      }

      FloatType z()     { return z_; }
      FloatType scale() { return scale_; }
      FloatType tp()    { return tp_; }
      FloatType hw()    { return hw_; }
      FloatType hn()    { return hn_; }
  };

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_STATISTICS_H
