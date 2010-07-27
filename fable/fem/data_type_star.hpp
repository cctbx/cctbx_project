#ifndef FEM_DATA_TYPE_STAR_HPP
#define FEM_DATA_TYPE_STAR_HPP

#include <fem/utils/int_types.hpp>

namespace fem {

  typedef bool logical_star_1;

  typedef utils::int8_t  integer_star_1;
  typedef utils::int16_t integer_star_2;
#if defined(_MSC_VER)
  typedef int            integer_star_4;
#else
  typedef utils::int32_t integer_star_4;
#endif
  typedef utils::int64_t integer_star_8;

  typedef float  real_star_4;
  typedef double real_star_8;

} // namespace fem

#endif // GUARD
