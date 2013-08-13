#ifndef SCITBX_BOOSTPYTHON_PYTHONITERATORRANGE_HPP_
#define SCITBX_BOOSTPYTHON_PYTHONITERATORRANGE_HPP_

#include <boost/python/object.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/class.hpp>

#include <boost/range/iterator_range.hpp>

namespace scitbx
{

namespace boost_python
{

template< class Iterator >
boost::python::object
as_iterator(Iterator const& begin, Iterator const& end)
{
  typedef boost::iterator_range< Iterator > range_type;
  return boost::python::iterator< range_type >()( range_type( begin, end ) );
}

template< class Iterator >
void
export_range_as_iterator(char const* name)
{
  typedef boost::iterator_range< Iterator > range_type;
  boost::python::class_< range_type >( name, boost::python::no_init )
    ;
}

} // namespace boost_python
} // namespace scitbx

#endif // SCITBX_BOOSTPYTHON_PYTHONITERATORRANGE_HPP_

