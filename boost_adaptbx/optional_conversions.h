#ifndef BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H
#define BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H

#include <boost/python/to_python_converter.hpp>
#include <boost/python/extract.hpp>
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

  template <typename T>
  struct from_python
  {
    from_python()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<boost::optional<T> >());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      if (obj_ptr != Py_None) {
        boost::python::extract<T> proxy(obj_ptr);
        if (!proxy.check()) return 0;
      }
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      boost::optional<T> value;
      if (obj_ptr != Py_None) {
        boost::python::extract<T> proxy(obj_ptr);
        value.reset(proxy());
      }
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<
          boost::optional<T> >*)
            data)->storage.bytes;
      new (storage) boost::optional<T>(value);
      data->convertible = storage;
    }
  };

  template <typename T>
  struct to_and_from_python
  {
    to_and_from_python()
    {
      to_python<T>();
      from_python<T>();
    }
  };

}} // namespace boost_adptbx::optional_conversions

#endif // BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H
