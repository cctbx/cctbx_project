/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Mar: Created, based on cctbx::maptbx::peak_histogram (rwgk)
 */

#ifndef SCITBX_HISTOGRAM_H
#define SCITBX_HISTOGRAM_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref_reductions.h>

namespace scitbx {

  //! Histogram of an array of scalar values.
  template <typename ValueType = double,
            typename CountType = long>
  class histogram
  {
    public:
      histogram() {}

      //! Computation of the histogram.
      template <typename DataType>
      histogram(
        af::const_ref<DataType> const& data,
        std::size_t n_slots=1000)
      :
        data_min_(0),
        data_max_(0),
        slot_width_(0),
        slots_(n_slots)
      {
        SCITBX_ASSERT(n_slots > 0);
        if (data.size() == 0) return;
        data_min_ = af::min(data);
        data_max_ = af::max(data);
        slot_width_ = (data_max_ - data_min_) / slots_.size();
        for(std::size_t i=0;i<data.size();i++) {
          ValueType d = data[i] - data_min_;
          std::size_t i_slot = 0;
          if (d != 0 && d >= slot_width_) {
                i_slot = std::size_t(d / slot_width_);
            if (i_slot >= slots_.size()) i_slot = slots_.size() - 1;
          }
          slots_[i_slot]++;
        }
      }

      //! Minimum value in data array passed to the constructor.
      ValueType
      data_min() const { return data_min_; }

      //! Maximum value in data array passed to the constructor.
      ValueType
      data_max() const { return data_max_; }

      //! Slot width used in the determination of the histogram.
      ValueType
      slot_width() const { return slot_width_; }

      //! Direct access to the array of counts.
      af::shared<CountType>
      slots() const { return slots_; }

      //! Determination of the cutoff value given a maximum number of points.
      ValueType
      get_cutoff(CountType const& max_points,
                 ValueType const& tolerance=1.e-4) const
      {
        CountType cum = 0;
        std::size_t i = slots_.size();
        for (; i; i--) {
          cum += slots_[i-1];
          if (cum > max_points) break;
        }
        return data_min_ + i * slot_width_ + slot_width_ * tolerance;
      }

    private:
      ValueType data_min_;
      ValueType data_max_;
      ValueType slot_width_;
      af::shared<CountType> slots_;
  };

} // namespace scitbx

#endif // SCITBX_HISTOGRAM_H
