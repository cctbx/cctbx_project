// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     May 2002: Created using code from phenix/translation.h (rwgk)
 */

#ifndef CCTBX_MAPS_GRIDDING_H
#define CCTBX_MAPS_GRIDDING_H

#include <cctbx/uctbx.h>
#include <cctbx/fftbx/gridding.h>

namespace cctbx { namespace maps {

  template <typename IndexType>
  IndexType
  determine_grid(const uctbx::UnitCell& UCell,
                 double max_Q,
                 double resolution_factor,
                 int max_prime,
                 const IndexType& mandatory_factors)
  {
    cctbx_assert(resolution_factor <= 0.5);
    double min_d = cctbx::uctbx::Q_as_d(max_Q);
    cctbx::af::int3 grid(
      UCell.MaxMillerIndices(min_d * 2 * resolution_factor));
    for(std::size_t i=0;i<3;i++) {
      grid[i] = 2 * grid[i] + 1;
    }
    return cctbx::fftbx::adjust_gridding_array(
      grid, max_prime, mandatory_factors);
  }

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_GRIDDING_H
