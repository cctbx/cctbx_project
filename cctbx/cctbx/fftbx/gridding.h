// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jan  2: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_FFTBX_GRIDDING_H
#define CCTBX_FFTBX_GRIDDING_H

#include <cctbx/fftbx/error.h>
#include <cctbx/fftbx/factorization.h>

namespace cctbx { namespace fftbx {

  template <typename IntegerType>
  bool
  check_max_prime(const IntegerType& max_prime, const IntegerType& N)
  {
    IntegerType RedN = N;
    detail::CountReduce(RedN, IntegerType(2));
    for (IntegerType factor = 3; RedN > 1; factor += 2) {
      if (factor > max_prime) return false;
      detail::CountReduce(RedN, factor);
    }
    return true;
  }

  template <typename IntegerType>
  IntegerType
  adjust_gridding(
    const IntegerType& min_grid,
    IntegerType max_prime,
    IntegerType mandatory_factor = 1)
  {
    if (max_prime < 2) max_prime = 0;
    if (mandatory_factor < 2) mandatory_factor = 1;
    IntegerType grid = (min_grid / mandatory_factor) * mandatory_factor;
    if (grid < min_grid) grid += mandatory_factor;
    if (max_prime) {
      if (!check_max_prime(max_prime, mandatory_factor)) {
        throw error(
          "adjust_gridding: mandatory_factor contains prime > max_prime");
      }
      while (!check_max_prime(max_prime, grid)) grid += mandatory_factor;
    }
    return grid;
  }

  template <typename IntegerArrayType>
  IntegerArrayType
  adjust_gridding_array(
    const IntegerArrayType& min_grid,
    const typename IntegerArrayType::value_type& max_prime)
  {
    IntegerArrayType result;
    for(std::size_t i=0;i<min_grid.size();i++) {
      result[i] = adjust_gridding(min_grid[i], max_prime);
    }
    return result;
  }

  template <typename IntegerArrayType>
  IntegerArrayType
  adjust_gridding_array(
    const IntegerArrayType& min_grid,
    const typename IntegerArrayType::value_type& max_prime,
    const IntegerArrayType& mandatory_factors)
  {
    IntegerArrayType result;
    for(std::size_t i=0;i<min_grid.size();i++) {
      result[i] = adjust_gridding(min_grid[i], max_prime,
                                  mandatory_factors[i]);
    }
    return result;
  }

}} // namespace cctbx::fftbx

#endif // CCTBX_FFTBX_GRIDDING_H
