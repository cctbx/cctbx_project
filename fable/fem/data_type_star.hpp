#ifndef FEM_DATA_TYPE_STAR_HPP
#define FEM_DATA_TYPE_STAR_HPP

#include <boost/cstdint.hpp>

namespace fem {

  typedef bool logical_star_1;

  typedef boost::int8_t integer_star_1;
  typedef boost::int16_t integer_star_2;
#if defined(_MSC_VER)
  typedef int            integer_star_4;
#else
  typedef boost::int32_t integer_star_4;
#endif
  typedef boost::int64_t integer_star_8;

  typedef float real_star_4;
  typedef double real_star_8;

} // namespace fem

#endif // GUARD
