#ifndef BOOST_ADAPTBX_CONVERTIBLE_H
#define BOOST_ADAPTBX_CONVERTIBLE_H

namespace boost_adaptbx {

template <class ConvertibleType, class ConvertedType,
          class ReportedConvertedType=ConvertedType>
struct convertible_type_to_python
{
  static PyObject *convert(ConvertibleType const &e) {
    using namespace boost::python;
    ConvertedType result = e;
    return incref(object(result).ptr());
  }

  static PyTypeObject const *get_pytype() {
    using namespace boost::python;
    return converter::registered<ReportedConvertedType>
              ::converters.to_python_target_type();
  }

  convertible_type_to_python() {
    using namespace boost::python;
    to_python_converter<
      ConvertibleType,
      convertible_type_to_python<
        ConvertibleType, ConvertedType, ReportedConvertedType>
      #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
        , true
      #endif
      >();
  }
};


}

#endif // GUARD
