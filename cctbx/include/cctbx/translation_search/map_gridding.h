/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified fragment from phenix/translation_search.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_MAP_GRIDDING_H
#define CCTBX_TRANSLATION_SEARCH_MAP_GRIDDING_H

#include <cctbx/maptbx/symmetry_flags.h>
#include <cctbx/maptbx/gridding.h>

namespace cctbx { namespace translation_search {

  template <typename GriddingTupleType = af::tiny<int, 3> >
  class map_gridding
  {
    public:
      map_gridding() {}

      map_gridding(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group_type const& sg_type,
        maptbx::symmetry_flags const& symmetry_flags,
        double resolution_factor,
        af::const_ref<miller::index<> > const& miller_indices_f_obs,
        int max_prime);

      const GriddingTupleType&
      target() const { return target_; }

      const GriddingTupleType&
      quarter() const { return quarter_; }

      const GriddingTupleType&
      eighth() const { return eighth_; }

    protected:
      GriddingTupleType target_;
      GriddingTupleType quarter_;
      GriddingTupleType eighth_;
  };

  template <typename GriddingTupleType>
  map_gridding<GriddingTupleType>
  ::map_gridding(
    uctbx::unit_cell const& unit_cell,
    sgtbx::space_group_type const& sg_type,
    maptbx::symmetry_flags const& symmetry_flags,
    double resolution_factor,
    af::const_ref<miller::index<> > const& miller_indices_f_obs,
    int max_prime)
  {
    GriddingTupleType
      sym_grid_factors = symmetry_flags.grid_factors(sg_type);
    double d_min = uctbx::d_star_sq_as_d(
      unit_cell.max_d_star_sq(miller_indices_f_obs));
    target_ = maptbx::determine_grid(
      unit_cell, d_min, resolution_factor, max_prime, sym_grid_factors);
    quarter_ = maptbx::determine_grid(
      unit_cell, d_min, 1./4, max_prime, target_);
    eighth_ = maptbx::determine_grid(
      unit_cell, d_min, 1./8, max_prime, target_);
  }

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_MAP_GRIDDING_H
