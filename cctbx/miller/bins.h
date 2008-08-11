#ifndef CCTBX_MILLER_BINS_H
#define CCTBX_MILLER_BINS_H

#include <cctbx/miller.h>
#include <cctbx/uctbx.h>
#include <cctbx/error.h>
#include <scitbx/math/linear_interpolation.h>

namespace cctbx { namespace miller {

  template <typename FloatType>
  inline
  FloatType
  sphere_volume(FloatType const& radius)
  {
    return scitbx::constants::four_pi * radius * radius * radius / 3;
  }

  class binning
  {
    public:
      binning() {}

      binning(
        uctbx::unit_cell const& unit_cell,
        std::size_t n_bins,
        double d_max,
        double d_min,
        double relative_tolerance=1.e-6)
      :
        unit_cell_(unit_cell)
      {
        init_limits(n_bins, d_max, d_min, relative_tolerance);
      }

      binning(
        uctbx::unit_cell const& unit_cell,
        std::size_t n_bins,
        af::const_ref<index<> > const& miller_indices,
        double d_max=0,
        double d_min=0,
        double relative_tolerance=1.e-6);

      binning(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<index<> > const& miller_indices,
        double d_max,
        double d_min,
        double d_star_sq_step);

      binning(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<double> const& limits)
      :
        unit_cell_(unit_cell),
        limits_(limits.begin(), limits.end())
      {}

      uctbx::unit_cell const&
      unit_cell() const { return unit_cell_; }

      std::size_t
      n_bins_used() const { return limits_.size() - 1; }

      std::size_t
      n_bins_all() const { return limits_.size() + 1; }

      std::size_t
      i_bin_d_too_large() const { return 0; }

      std::size_t
      i_bin_d_too_small() const { return limits_.size(); }

      double
      d_max() const { return bin_d_min(1); }

      double
      d_min() const { return bin_d_min(n_bins_all()-1); }

      af::double2
      bin_d_range(std::size_t i_bin) const;

      double
      bin_d_min(std::size_t i) const;

      af::shared<double> const&
      limits() const { return limits_; }

      std::size_t
      get_i_bin(double d_star_sq) const;

      std::size_t
      get_i_bin(index<> const& h) const
      {
        return get_i_bin(unit_cell_.d_star_sq(h));
      }

    protected:
      void
      init_limits(
        std::size_t n_bins,
        double d_star_sq_min,
        double d_star_sq_max,
        double relative_tolerance);

      void
      init_limits_d_star_sq_step(
        double d_min,
        double d_max,
        double d_star_sq_step);


      uctbx::unit_cell unit_cell_;
      af::shared<double> limits_;
  };

  class binner : public binning
  {
    public:
      binner() {}

      binner(
        binning const& bng,
        af::shared<index<> > const& miller_indices);

      af::shared<index<> > const&
      miller_indices() const { return miller_indices_; }

      af::shared<std::size_t> const&
      bin_indices() const { return bin_indices_; }

      std::size_t
      count(std::size_t i_bin) const;

      af::shared<std::size_t>
      counts() const;

      af::shared<bool>
      selection(std::size_t i_bin) const;

      af::shared<std::size_t>
      array_indices(std::size_t i_bin) const;

      template <typename FloatType>
      af::shared<FloatType>
      bin_centers(FloatType const& d_star_power) const;

      template <typename FloatType>
      af::shared<FloatType>
      interpolate(
        af::const_ref<FloatType> const& values,
        FloatType const& d_star_power) const;

    private:
      af::shared<index<> > miller_indices_;
      af::shared<std::size_t> bin_indices_;
  };

  namespace detail {

    template <typename FloatType>
    inline
    FloatType
    d_star_to_the(
      FloatType const& d_star_sq,
      FloatType const& d_star_power)
    {
      if (d_star_power == 1) return std::sqrt(d_star_sq);
      if (d_star_power == 2) return d_star_sq;
      return std::pow(d_star_sq, d_star_power/2);
    }

  } // namespace detail

  template <typename FloatType>
  af::shared<FloatType>
  binner::bin_centers(FloatType const& d_star_power) const
  {
    std::size_t l_bin = this->i_bin_d_too_large() + 1;
    std::size_t s_bin = this->i_bin_d_too_small() - 1;
    af::shared<FloatType> result((af::reserve(this->n_bins_used())));
    FloatType a = detail::d_star_to_the(this->limits_[0], d_star_power);
    for(std::size_t i_bin=l_bin;i_bin<=s_bin;i_bin++) {
      FloatType b = detail::d_star_to_the(this->limits_[i_bin], d_star_power);
      result.push_back(a+(b-a)/2);
      a = b;
    }
    return result;
  }

  template <typename FloatType>
  af::shared<FloatType>
  binner::interpolate(
    af::const_ref<FloatType> const& values,
    FloatType const& d_star_power) const
  {
    af::const_ref<index<> > miller_indices = miller_indices_.const_ref();
    CCTBX_ASSERT(miller_indices.size() == bin_indices_.size());
    CCTBX_ASSERT(values.size() == this->n_bins_used());
    af::shared<FloatType> result((af::reserve(miller_indices.size())));
    if (d_star_power == 0 || values.size() == 1) {
      for(std::size_t i=0;i<miller_indices.size();i++) {
        result.push_back(values[bin_indices_[i]-1]);
      }
    }
    else {
      af::shared<FloatType> bin_ctrs_mem = bin_centers(d_star_power);
      af::const_ref<FloatType> bin_ctrs = bin_ctrs_mem.const_ref();
      for(std::size_t i=0;i<miller_indices.size();i++) {
        std::size_t i_bin = bin_indices_[i];
        if (   i_bin == this->i_bin_d_too_large()
            || i_bin == this->i_bin_d_too_small()) {
          throw error("Miller index outside binned range.");
        }
        FloatType x = detail::d_star_to_the(
          this->unit_cell_.d_star_sq(miller_indices[i]),
          d_star_power);
        FloatType c = bin_ctrs[i_bin-1];
        std::size_t l_bin;
        if (x < c) l_bin = std::max(i_bin-1, this->i_bin_d_too_large()+1);
        else l_bin = i_bin;
        std::size_t s_bin = std::min(l_bin+1, this->i_bin_d_too_small()-1);
        l_bin = s_bin - 1;
        CCTBX_ASSERT(l_bin > this->i_bin_d_too_large());
        result.push_back(scitbx::math::linear_interpolation(
          x,
          bin_ctrs[l_bin-1], bin_ctrs[s_bin-1],
          values[l_bin-1], values[s_bin-1]));
      }
    }
    return result;
  }






}} // namespace cctbx::miller

#endif // CCTBX_MILLER_BINS_H
