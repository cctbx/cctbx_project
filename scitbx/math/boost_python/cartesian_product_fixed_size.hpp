#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/stl_iterator.hpp>

#include <scitbx/math/cartesian_product_fixed_size.hpp>

namespace scitbx
{

namespace math
{

namespace cartesian_product
{

namespace python
{

template< typename Iterator >
struct iterated_range_wrappers
{
  static void wrap(const char* name)
  {
    using namespace boost::python;

    typedef iterated_range< Iterator > range_type;

    class_< range_type >( name, no_init )
      .def(
        "__iter__",
        boost::python::range( &range_type::begin, &range_type::end )
        )
      .def( "__len__", &range_type::length )
      ;
  }
};

} // namespace python
} // namespace cartesian_product
} // namespace math
} // namespace scitbx

