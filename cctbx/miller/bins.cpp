/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Refactored parts of miller/miller_lib.cpp (rwgk)
 */

#include <cctbx/miller/bins.h>
#include <cctbx/error.h>

namespace cctbx { namespace miller {

  binning::binning(
    uctbx::unit_cell const& unit_cell,
    std::size_t n_bins,
    af::const_ref<index<> > const& miller_indices,
    double d_max,
    double d_min,
    double relative_tolerance)
  :
    unit_cell_(unit_cell)
  {
    if (!(d_max || d_min)) {
      af::double2 min_max_q = unit_cell.min_max_d_star_sq(miller_indices);
      if (!d_max && min_max_q[0]) d_max = 1 / std::sqrt(min_max_q[0]);
      if (!d_min && min_max_q[1]) d_min = 1 / std::sqrt(min_max_q[1]);
    }
    init_limits(n_bins, d_max, d_min, relative_tolerance);
  }

  void
  binning::init_limits(
    std::size_t n_bins,
    double d_max,
    double d_min,
    double relative_tolerance)
  {
    CCTBX_ASSERT(n_bins > 0);
    CCTBX_ASSERT(d_max >= 0);
    CCTBX_ASSERT(d_min > 0);
    CCTBX_ASSERT(d_max == 0 || d_min < d_max);
    double d_star_sq_min = 0;
    if (d_max) d_star_sq_min = 1 / (d_max * d_max);
    double     d_star_sq_max = 1 / (d_min * d_min);
    double span = d_star_sq_max - d_star_sq_min;
    d_star_sq_max += span * relative_tolerance;
    d_star_sq_min -= span * relative_tolerance;
    if (d_star_sq_min < 0) d_star_sq_min = 0;
    double r_low = std::sqrt(d_star_sq_min);
    double r_high = std::sqrt(d_star_sq_max);
    double volume_low = sphere_volume(r_low);
    double volume_per_bin = (sphere_volume(r_high) - volume_low) / n_bins;
    limits_.push_back(d_star_sq_min);
    for(std::size_t i_bin=1;i_bin<n_bins;i_bin++) {
      double r_sq_i = std::pow(
        (volume_low + i_bin * volume_per_bin) * 3 / scitbx::constants::four_pi,
        2/3.);
      limits_.push_back(r_sq_i);
    }
    limits_.push_back(d_star_sq_max);
  }

  af::double2
  binning::bin_d_range(std::size_t i_bin) const
  {
    return af::double2(bin_d_min(i_bin), bin_d_min(i_bin+1));
  }

  double
  binning::bin_d_min(std::size_t i_bin) const
  {
    if (i_bin == 0) return -1;
    if (i_bin == n_bins_all()) return -1;
    if (i_bin > n_bins_all()) throw error_index();
    return uctbx::d_star_sq_as_d(limits_[i_bin - 1]);
  }

  std::size_t
  binning::get_i_bin(double d_star_sq) const
  {
    if (d_star_sq < limits_[0]) return 0;
    std::size_t i = 1;
    for(;i<limits_.size();i++) {
      if (d_star_sq < limits_[i]) return i;
    }
    return i;
  }

  binner::binner(
    binning const& bng,
    af::shared<index<> > const& miller_indices)
  :
    binning(bng),
    miller_indices_(miller_indices)
  {
    af::const_ref<index<> > mi = miller_indices_.const_ref();
    bin_indices_.reserve(mi.size());
    for(std::size_t i=0;i<mi.size();i++) {
      bin_indices_.push_back(this->get_i_bin(mi[i]));
    }
  }

  std::size_t
  binner::count(std::size_t i_bin) const
  {
    CCTBX_ASSERT(i_bin < this->n_bins_all());
    std::size_t result = 0;
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      if (bin_indices_[i] == i_bin) result++;
    }
    return result;
  }

  af::shared<std::size_t>
  binner::counts() const
  {
    af::shared<std::size_t> result(this->n_bins_all());
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      std::size_t i_bin = bin_indices_[i];
      CCTBX_ASSERT(i_bin < result.size());
      result[i_bin]++;
    }
    return result;
  }

  af::shared<bool>
  binner::selection(std::size_t i_bin) const
  {
    CCTBX_ASSERT(i_bin < this->n_bins_all());
    af::shared<bool> flags((af::reserve(bin_indices_.size())));
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      flags.push_back(bin_indices_[i] == i_bin);
    }
    return flags;
  }

  af::shared<std::size_t>
  binner::array_indices(std::size_t i_bin) const
  {
    CCTBX_ASSERT(i_bin < this->n_bins_all());
    af::shared<std::size_t> result((af::reserve(bin_indices_.size())));
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      if (bin_indices_[i] == i_bin) result.push_back(i);
    }
    return result;
  }

}} // namespace cctbx::miller
