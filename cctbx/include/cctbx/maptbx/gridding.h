/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_GRIDDING_H
#define CCTBX_MAPTBX_GRIDDING_H

#include <cctbx/maptbx/symmetry_flags.h>
#include <scitbx/fftpack/gridding.h>
#include <scitbx/array_family/tiny_algebra.h>

namespace cctbx { namespace maptbx {

  template <typename IndexValueType>
  af::tiny<IndexValueType, 3>
  determine_gridding(
    uctbx::unit_cell const& unit_cell,
    double d_min,
    double resolution_factor,
    af::tiny<IndexValueType, 3> const& mandatory_factors,
    IndexValueType max_prime=5,
    bool assert_shannon_sampling=true)
  {
    CCTBX_ASSERT(d_min > 0);
    if (assert_shannon_sampling) {
      CCTBX_ASSERT(resolution_factor <= 0.5);
    }
    af::tiny<IndexValueType, 3>
      grid(unit_cell.max_miller_indices(d_min * 2 * resolution_factor));
    grid *= IndexValueType(2);
    grid += IndexValueType(1);
    return scitbx::fftpack::adjust_gridding_array(
      grid, max_prime, mandatory_factors);
  }

  // Not available in Python.
  template <typename IndexValueType>
  af::tiny<IndexValueType, 3>
  determine_gridding(
    uctbx::unit_cell const& unit_cell,
    double d_min,
    double resolution_factor,
    IndexValueType max_prime=5,
    bool assert_shannon_sampling=true)
  {
    return determine_gridding(
      unit_cell, d_min, resolution_factor,
      af::tiny<IndexValueType, 3>(1,1,1),
      max_prime, assert_shannon_sampling);
  }

  template <typename IndexValueType>
  af::tiny<IndexValueType, 3>
  determine_gridding(
    uctbx::unit_cell const& unit_cell,
    double d_min,
    double resolution_factor,
    maptbx::symmetry_flags const& symmetry_flags,
    sgtbx::space_group_type const& space_group_type,
    IndexValueType max_prime=5,
    bool assert_shannon_sampling=true)
  {
    typedef IndexValueType i_v_t;
    sgtbx::structure_seminvariant ss(space_group_type.group());
    af::tiny<i_v_t, 3> mandatory_factors(1,1,1);
    if (symmetry_flags.use_structure_seminvariants()) {
      mandatory_factors = ss.refine_gridding(mandatory_factors);
    }
    sgtbx::space_group sub_space_group = symmetry_flags.select_sub_space_group(
      space_group_type);
    mandatory_factors = sub_space_group.refine_gridding(mandatory_factors);
    af::tiny<i_v_t, 3>
      grid = determine_gridding(
        unit_cell, d_min, resolution_factor,
        mandatory_factors, max_prime, assert_shannon_sampling);
    std::size_t best_size = 0;
    af::tiny<i_v_t, 3> best_grid(0,0,0);
    i_v_t g_limit = af::max(grid) + 1;
    af::tiny<i_v_t, 3> loop;
    for(loop[0]=grid[0];loop[0]<g_limit;loop[0]+=mandatory_factors[0])
    for(loop[1]=grid[1];loop[1]<g_limit;loop[1]+=mandatory_factors[1])
    for(loop[2]=grid[2];loop[2]<g_limit;loop[2]+=mandatory_factors[2]) {
      af::tiny<i_v_t, 3> trial_grid = scitbx::fftpack::adjust_gridding_array(
        loop, max_prime, mandatory_factors);
      if (symmetry_flags.use_structure_seminvariants()) {
        trial_grid = ss.refine_gridding(trial_grid);
      }
      trial_grid = sub_space_group.refine_gridding(trial_grid);
      CCTBX_ASSERT(scitbx::fftpack::adjust_gridding_array(
        trial_grid, max_prime, mandatory_factors).all_eq(trial_grid));
      if (best_size == 0 && trial_grid.all_eq(grid)) {
        return grid;
      }
      std::size_t trial_size = 1;
      for(std::size_t i=0;i<3;i++) trial_size *= trial_grid[i];
      CCTBX_ASSERT(trial_size > 0);
      if (best_size == 0 || trial_size < best_size) {
        best_grid = trial_grid;
        best_size = trial_size;
      }
    }
    return best_grid;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_GRIDDING_H
