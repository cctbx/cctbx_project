/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: modified copy of cctbx/math/array_utils.h (rwgk)
     2002 Mar: modified copy of cctbx/maps/utils.h (rwgk)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_STATISTICS_H
#define CCTBX_MAPTBX_STATISTICS_H

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/loops.h>
#include <cctbx/math/utils.h>
#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace maptbx {

  //! Determines simple map statistics.
  template <typename FloatType = double>
  class statistics
  {
    public:
      //! Default constructor. Data members are not initialized!
      statistics() {}

      //! Computes the statistics.
      template <typename OtherFloatType>
      statistics(af::const_ref<OtherFloatType, af::flex_grid<> > const& map)
      {
        CCTBX_ASSERT(map.accessor().focus_size_1d() > 0);
        if (!map.accessor().is_padded()) {
          min_ = af::min(map);
          max_ = af::max(map);
          mean_ = af::mean(map);
          mean_sq_ = af::mean_sq(map);
        }
        else {
          min_ = max_ = map[0];
          mean_ = mean_sq_ = FloatType(0);
          typedef typename af::flex_grid<>::index_type index_type;
          af::flex_grid<> zero_based = map.accessor().shift_origin();
          af::nested_loop<index_type> loop(zero_based.focus());
          std::size_t n = 0;
          for (index_type const& pt = loop(); !loop.over(); loop.incr()) {
            FloatType v = map[zero_based(pt)];
            math::update_min(min_, v);
            math::update_max(max_, v);
            mean_ += v;
            mean_sq_ += v * v;
            n++;
          }
          mean_ /= FloatType(n);
          mean_sq_ /= FloatType(n);
        }
        sigma_ = mean_sq_ - mean_ * mean_;
        if (sigma_ < FloatType(0)) sigma_ = 0;
        sigma_ = std::sqrt(sigma_);
      }

      FloatType
      min() const { return min_; }

      FloatType
      max() const { return max_; }

      FloatType
      mean() const { return mean_; }

      FloatType
      mean_sq() const { return mean_sq_; }

      FloatType
      sigma() const { return sigma_; }

    protected:
      FloatType min_;
      FloatType max_;
      FloatType mean_;
      FloatType mean_sq_;
      FloatType sigma_;
  };

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_STATISTICS_H
