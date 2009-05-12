#ifndef SCITBX_OPTIONAL_COPY_CONVERSION_H
#define SCITBX_OPTIONAL_COPY_CONVERSION_H

#include <boost/python/to_python_converter.hpp>
#include <boost/python/extract.hpp>

#include <scitbx/optional_copy.h>

namespace scitbx { namespace boost_python {

namespace optional_copy_conversions {

  template <class ValueType>
  struct to_python
  {
    static PyObject *convert(optional_copy<ValueType> const &value) {
      using namespace boost::python;
      object result;
      if (value.get()) result = object(value.get());
      return incref(result.ptr());
    }

    static PyTypeObject const *get_pytype() {
      return boost::python::converter::registered<ValueType>::converters
        .to_python_target_type();
    }

    to_python() {
      boost::python::to_python_converter<
        optional_copy<ValueType>,
        to_python
        #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
        #endif
        >();
    }
  };

}}} // scitbx::boost_python::optional_copy_conversions

#endif // GUARD
