/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_BINS_H
#define CCTBX_MILLER_BINS_H

#include <cctbx/miller.h>
#include <cctbx/uctbx.h>

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

      uctbx::unit_cell unit_cell_;
      af::shared<double> limits_;
  };

  class binner : public binning
  {
    public:
      binner() {}

      binner(
        binning const& bng,
        af::const_ref<index<> > const& miller_indices);

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

    private:
      af::shared<std::size_t> bin_indices_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_BINS_H
