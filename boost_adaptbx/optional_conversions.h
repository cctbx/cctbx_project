#ifndef BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H
#define BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H

#include <boost/python/to_python_converter.hpp>
#include <boost/python/object.hpp>
#include <boost/optional.hpp>

namespace boost_adaptbx { namespace optional_conversions {

  namespace detail {

    template <typename T>
    struct to_python
    {
      static PyObject*
      convert(
        boost::optional<T> const& value)
      {
        if (value) {
          return boost::python::incref(boost::python::object(*value).ptr());
        }
        return boost::python::incref(Py_None);
      }
    };

  } // namespace detail

  template <typename T>
  struct to_python
  {
    to_python()
    {
      boost::python::to_python_converter<
        boost::optional<T>,
        detail::to_python<T> >();
    }
  };

}} // namespace boost_adptbx::optional_conversions

#endif // BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H
