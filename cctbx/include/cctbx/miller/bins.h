// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_BINS_H
#define CCTBX_MILLER_BINS_H

#include <cctbx/uctbx.h>
#include <cctbx/miller.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace miller {

  template <typename FloatType>
  inline
  FloatType
  sphere_volume(FloatType const& radius)
  {
    return constants::four_pi * radius * radius * radius / 3;
  }

  class binning
  {
    public:
      binning() {}

      binning(
        uctbx::UnitCell const& unit_cell,
        af::shared<Index> miller_indices,
        std::size_t n_bins,
        double d_max = 0,
        double d_min = 0);

      binning(
        uctbx::UnitCell const& unit_cell,
        std::size_t n_bins,
        double d_max,
        double d_min);

      uctbx::UnitCell const& unit_cell() const { return unit_cell_; }

      std::size_t n_bins() const { return limits_.size() - 1; }

      double d(std::size_t i) const
      {
        if (i >= limits_.size()) throw error_index();
        return 1 / std::sqrt(limits_[i]);
      }

      double d_min() const { return d(n_bins()); }

      af::shared<double> limits() const { return limits_; }

      std::size_t
      get_bin_index(double d_star_sq, double relative_tolerance=1.e-6) const;

      std::size_t
      get_bin_index(Index const& h, double relative_tolerance=1.e-6) const
      {
        return get_bin_index(unit_cell_.Q(h), relative_tolerance);
      }

    protected:
      void init_limits(
        double d_star_sq_min, double d_star_sq_max, std::size_t n_bins);

      uctbx::UnitCell unit_cell_;
      af::shared<double> limits_;
      double span_;
  };

  class binner
  {
    public:
      binner() {}

      binner(binning const& bng, af::shared<Index> miller_indices);

      af::shared<std::size_t> bin_indices() const { return bin_indices_; }

      af::shared<bool> bin_selection(std::size_t i_bin) const;

    private:
      af::shared<std::size_t> bin_indices_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_BINS_H
