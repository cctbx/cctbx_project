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
      get_j_bin(double d_star_sq, double relative_tolerance=1.-6) const
      {
        double abs_tolerance = relative_tolerance * span_;
        if (d_star_sq > limits_[0] + abs_tolerance) return 0;
        std::size_t i = 1;
        for(;i<limits_.size();i++) {
          if (d_star_sq > limits_[i]) return i;
        }
        if (d_star_sq > limits_[i] + abs_tolerance) i--;
        return i;
      }
      std::size_t
      get_j_bin(Index const& h, double relative_tolerance=1.-6) const
      {
        return get_j_bin(unit_cell_.Q(h), relative_tolerance);
      }

    protected:
      void init_limits(
        double d_star_sq_min, double d_star_sq_max, std::size_t n_bins);

      uctbx::UnitCell unit_cell_;
      af::shared<double> limits_;
      double span_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_BINS_H
