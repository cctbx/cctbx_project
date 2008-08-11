#ifndef SCITBX_FFTPACK_GRIDDING_H
#define SCITBX_FFTPACK_GRIDDING_H

#include <scitbx/error.h>
#include <scitbx/fftpack/factorization.h>

namespace scitbx { namespace fftpack {

  template <typename IntegerType>
  bool
  check_max_prime(const IntegerType& max_prime, const IntegerType& n)
  {
    IntegerType red_n = n;
    detail::count_reduce(red_n, IntegerType(2));
    for (IntegerType factor = 3; red_n > 1; factor += 2) {
      if (factor > max_prime) return false;
      detail::count_reduce(red_n, factor);
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

  // XXX this function should eventually be obsolete
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

  // XXX this function should eventually be obsolete
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

  template <typename IntegerArrayType>
  IntegerArrayType
  adjust_gridding_array_flex(
    IntegerArrayType const& min_grid,
    typename IntegerArrayType::value_type const& max_prime)
  {
    IntegerArrayType result;
    for(std::size_t i=0;i<min_grid.size();i++) {
      result.push_back(adjust_gridding(min_grid[i], max_prime));
    }
    return result;
  }

  template <typename IntegerArrayType>
  IntegerArrayType
  adjust_gridding_array_flex(
    IntegerArrayType const& min_grid,
    typename IntegerArrayType::value_type const& max_prime,
    IntegerArrayType const& mandatory_factors)
  {
    if (min_grid.size() != mandatory_factors.size()) {
      throw error(
        "adjust_gridding_array: min_grid.size() != mandatory_factors.size()");
    }
    IntegerArrayType result;
    for(std::size_t i=0;i<min_grid.size();i++) {
      result.push_back(adjust_gridding(min_grid[i], max_prime,
                                       mandatory_factors[i]));
    }
    return result;
  }

}} // namespace scitbx::fftpack

#endif // SCITBX_FFTPACK_GRIDDING_H
