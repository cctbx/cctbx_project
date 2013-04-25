#ifndef BOOST_ADAPTBX_BOOST_RANGE_PYTHON_H
#define BOOST_ADAPTBX_BOOST_RANGE_PYTHON_H

#include <boost/range.hpp>

#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>

namespace boost_adaptbx
{

namespace python
{

template< typename Range >
struct generic_range_wrapper
{
  static void wrap(const char* name)
  {
    using namespace boost::python;

    class_< Range >( name, no_init )
      .def(
        "__iter__",
        boost::python::iterator< Range >()
        )
      .def( "__len__", boost::distance< Range > )
      .def( "empty", boost::empty< Range > )
      ;
  }
};

} // namespace python
} // namespace boost_adaptbx

#endif

