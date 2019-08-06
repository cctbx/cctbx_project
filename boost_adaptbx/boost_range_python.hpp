#ifndef BOOST_ADAPTBX_BOOST_RANGE_PYTHON_H
#define BOOST_ADAPTBX_BOOST_RANGE_PYTHON_H

#include <boost/range.hpp>

#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/converter/registry.hpp>

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

    boost::python::type_info info = boost::python::type_id< Range >();
    const converter::registration* reg = converter::registry::query( info );

    if ( reg == NULL || reg->m_to_python == NULL )
    {
      class_< Range >( name, no_init )
        .def(
          "__iter__",
          boost::python::iterator< Range >()
          )
        .def( "__len__", boost::size< Range > )
        .def( "empty", boost::empty< Range > )
        ;
    }
  }
};

} // namespace python
} // namespace boost_adaptbx

#endif

