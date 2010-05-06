#ifndef BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H
#define BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H

#include <boost/python/to_python_converter.hpp>
#include <boost/python/extract.hpp>

namespace boost_adaptbx { namespace optional_conversions {

  template <typename OptionalType>
  struct to_python
  {
    typedef typename OptionalType::value_type value_type;

    static PyObject*
    convert(
      OptionalType const& value)
    {
      if (value) {
        return boost::python::incref(boost::python::object(*value).ptr());
      }
      return boost::python::incref(Py_None);
    }

    static const PyTypeObject* get_pytype()
    {
      return boost::python::converter::registered<value_type>::converters
        .to_python_target_type();
    }

    to_python()
    {
      boost::python::to_python_converter<
        OptionalType,
        to_python<OptionalType>
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
#endif
        >();
    }
  };

  template <typename OptionalType>
  struct from_python
  {
    typedef typename OptionalType::value_type value_type;

    from_python()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<OptionalType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
      , &boost::python::converter::expected_pytype_for_arg<
          value_type>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
      if (obj_ptr != Py_None) {
        boost::python::extract<value_type> proxy(obj_ptr);
        if (!proxy.check()) return 0;
      }
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      OptionalType value;
      if (obj_ptr != Py_None) {
        boost::python::extract<value_type> proxy(obj_ptr);
        value = proxy();
      }
      void* storage = (
        (boost::python::converter::rvalue_from_python_storage<
          OptionalType>*)
            data)->storage.bytes;
      new (storage) OptionalType(value);
      data->convertible = storage;
    }
  };

  template <typename OptionalType>
  struct to_and_from_python
  {
    to_and_from_python()
    {
      to_python<OptionalType>();
      from_python<OptionalType>();
    }
  };

}} // namespace boost_adaptbx::optional_conversions

#endif // BOOST_ADAPTBX_OPTIONAL_CONVERSIONS_H
