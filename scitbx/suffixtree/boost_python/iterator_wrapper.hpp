#ifndef SUFFIXTREE_PYTHON_ITERATOR_WRAPPER_HPP_
#define SUFFIXTREE_PYTHON_ITERATOR_WRAPPER_HPP_

#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/def.hpp>

namespace scitbx
{
namespace suffixtree
{
namespace python
{

boost::python::object passthrough(boost::python::object const& o )
{
    return o;
}

template<
  typename InputIterator
  >
struct python_iterator 
{
  typedef InputIterator iterator;
  typedef typename InputIterator::value_type value_type;

  iterator pos_;
  iterator end_;

  python_iterator(iterator const& begin, iterator const& end)
  : pos_( begin ), end_( end )
  {}

  value_type next()
  {
      if ( pos_ == end_ )
      {
        PyErr_SetString( PyExc_StopIteration, "" );
        boost::python::throw_error_already_set();
      }

      return *( pos_++ );
  }

  static void wrap(const char* name)
  {
      boost::python::class_< python_iterator >( name, boost::python::no_init )
        .def( "next", &python_iterator::next )
        .def( "__iter__", passthrough )
        ;
  }
};

} // namespace python
} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_PYTHON_ITERATOR_WRAPPER_HPP_
