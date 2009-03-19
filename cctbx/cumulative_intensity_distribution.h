#ifndef CCTBX_CUMULATIVE_INTENSITY_DISTRIBUTION_H
#define CCTBX_CUMULATIVE_INTENSITY_DISTRIBUTION_H

#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>

namespace cctbx {

  template <typename FloatType>
  struct bin
  {
    // z = I/<I>
    bin(FloatType z_min_, FloatType z_max_) {
      count = 0;
      z_min = z_min_;
      z_max = z_max_;
    }
    int count;
    FloatType z_min;
    FloatType z_max;
  };

  template <typename FloatType>
  class cumulative_intensity
  // As described by  Howells, Phillips and Rogers, Acta Cryst. (1950). 3, 210
  {
  public:
    //! Constructor
    cumulative_intensity(
      af::const_ref<FloatType> const &f_sq,
      af::const_ref<FloatType> const &d_spacings,
      af::const_ref<FloatType> const &mean_f_sq_,
      af::const_ref<FloatType> const &bin_d_max_,
      af::shared<miller::index<> > const &indices)
    :
      mean_f_sq(mean_f_sq_),
      bin_d_max(bin_d_max_)
    {
      CCTBX_ASSERT(f_sq.size() == d_spacings.size());
      CCTBX_ASSERT(f_sq.size() == indices.size());
      CCTBX_ASSERT(mean_f_sq.size() == bin_d_max.size());

      const int n_bins = get_bin_count();
      af::shared<bin<FloatType> > binner;
      // setup binner
      for(int i=0;i<n_bins;) {
        FloatType z_min = i/FloatType(n_bins);
        FloatType z_max = ++i/FloatType(n_bins);
        binner.push_back(bin<FloatType>(z_min, z_max));
      }
      // loop over all reflections
      for(std::size_t i=0;i<indices.size();i++) {
        miller::index<> const &h = indices[i];
        FloatType i_over_mean_i = f_sq[i]/get_mean_f_sq(d_spacings[i]);
        // loop over bins
        for(std::size_t j=0;j<binner.size();j++) {
          FloatType const &z = i_over_mean_i;
          if (z < binner[j].z_max && z > binner[j].z_min) {
            binner[j].count ++;
            break;
          }
        }
      }
      // loop over bins again to get cumulative intensities
      FloatType cumulative_total = 0.0;
      for(std::size_t i=0;i<binner.size();i++) {
        x_.push_back(binner[i].z_max);
        cumulative_total += binner[i].count;
        y_.push_back(cumulative_total/FloatType(f_sq.size()));
      }
    }

    inline const af::shared<FloatType>& x() const { return x_; }
    inline const af::shared<FloatType>& y() const { return y_; }

  private:
    inline int get_bin_count() const { return mean_f_sq.size(); }

    FloatType get_mean_f_sq(FloatType const &d_spacing) const {
      for(std::size_t i=0;i<get_bin_count();i++) {
        if (d_spacing >= bin_d_max[i])
          return mean_f_sq[i];
      }
      throw std::runtime_error(
        "Unexpected d spacing, no bin found");
    }

    // result arrays
    af::shared<FloatType> x_;
    af::shared<FloatType> y_;

    af::const_ref<FloatType> mean_f_sq;
    af::const_ref<FloatType> bin_d_max;
  };
} // cctbx

#endif
