/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Moved from cctbx/maps to cctbx/maptbx (rwgk)
     2002 May: Created using code from phenix/translation.h (rwgk)
 */

#ifndef CCTBX_MAPTBX_GRIDDING_H
#define CCTBX_MAPTBX_GRIDDING_H

#include <cctbx/uctbx.h>
#include <scitbx/fftpack/gridding.h>
#include <cctbx/error.h>

namespace cctbx { namespace maptbx {

  template <typename IndexType>
  IndexType
  determine_grid(uctbx::unit_cell const& unit_cell,
                 double d_min,
                 double resolution_factor,
                 int max_prime,
                 IndexType const& mandatory_factors)
  {
    CCTBX_ASSERT(mandatory_factors.size() == 3);
    CCTBX_ASSERT(resolution_factor <= 0.5);
    IndexType grid(af::adapt(
      unit_cell.max_miller_indices(d_min * 2 * resolution_factor)));
    for(std::size_t i=0;i<grid.size();i++) {
      grid[i] = 2 * grid[i] + 1;
    }
    return scitbx::fftpack::adjust_gridding_array(
             grid, max_prime, mandatory_factors);
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_GRIDDING_H
