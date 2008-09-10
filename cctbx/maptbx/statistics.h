#ifndef CCTBX_MAPTBX_STATISTICS_H
#define CCTBX_MAPTBX_STATISTICS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/accumulators.h>
#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace maptbx {

  namespace details {
    using namespace scitbx::math::accumulator;

    template <typename FloatType = double>
    struct statistics_traits {
      typedef min_max_accumulator<FloatType,
                mean_variance_accumulator<FloatType,
                  enumerated_accumulator<FloatType> > >
              accumulator_t;
    };
  }

  //! Determines simple map statistics.
  template <typename FloatType = double>
  class statistics
  {
    public:
      //! Default constructor. Data members are not initialized!
      statistics() : stats(0) {}

      //! Computes the statistics.
      template <typename OtherFloatType>
      statistics(af::const_ref<OtherFloatType, af::flex_grid<> > const& map)
        : stats(0)
      {
        CCTBX_ASSERT(map.accessor().focus_size_1d() > 0);
        if (!map.accessor().is_padded()) {
          stats = accumulator_t(map[0]);
          for(std::size_t i=1; i < map.size(); ++i) stats(map[i]);
        }
        else {
          typedef typename af::flex_grid<>::index_type index_type;
          af::flex_grid<> zero_based = map.accessor().shift_origin();
          af::nested_loop<index_type> iter(zero_based.focus());
          stats = accumulator_t(map[zero_based(iter())]);
          while (iter.incr()) stats(map[zero_based(iter())]);
        }
      }

      FloatType
      min() const { return stats.min(); }

      FloatType
      max() const { return stats.max(); }

      FloatType
      mean() const { return stats.mean(); }

      FloatType
      mean_sq() const { return stats.mean_squares(); }

      FloatType
      sigma() const { return stats.biased_standard_deviation(); }

    protected:
      typedef typename details::statistics_traits<FloatType>::accumulator_t
              accumulator_t;
      accumulator_t stats;
  };

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_STATISTICS_H
