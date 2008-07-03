#ifndef BOOST_ADAPTBX_TUPLE_CONVERSION_H
#define BOOST_ADAPTBX_TUPLE_CONVERSION_H

#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/refcount.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/tuple/tuple.hpp>

/* From a recipe given by David Abraham on C++ sig
( http://mail.python.org/pipermail/c++-sig/2006-January/010067.html )
*/

namespace boost_adaptbx { namespace tuple_conversion {

namespace detail {

template <class T>
struct to_python
{
  template <class Head, class Tail>
  static
  inline
  boost::python::tuple
  tuple_to_python(boost::tuples::cons<Head, Tail> const& x) {
      boost::python::tuple head = boost::python::make_tuple(x.get_head());
      boost::python::tuple tail = tuple_to_python(x.get_tail());
      return boost::python::tuple(head + tail);
  }

  static
  inline
  boost::python::tuple
  tuple_to_python(boost::tuples::null_type) {
      return boost::python::tuple();
  }

  static PyObject *
  convert(T const& x) {
      return boost::python::incref(tuple_to_python(x).ptr());
  }

  static const PyTypeObject *
  get_pytype() {
    return boost::python::converter::registered<T>::converters
            .to_python_target_type();
  }
};

} // namespace detail

template<class T>
struct to_python
{
  to_python() {
    boost::python::to_python_converter<T, detail::to_python<T>
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
#endif
        >();
  }
};


}}

#endif // GUARD
