#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost_adaptbx/optional_conversions.h>
#include <boost_adaptbx/type_id_eq.h>
#include <boost/optional.hpp>

namespace {

  boost::optional<std::size_t>
  exercise(
    boost::optional<double> const& value)
  {
    if (value) {
      if (*value == 13) return boost::optional<std::size_t>();
      return boost::optional<std::size_t>(
        static_cast<std::size_t>((*value)*3));
    }
    return boost::optional<std::size_t>(42);
  }

#if !defined(BOOST_NO_STD_WSTRING) && defined(Py_USING_UNICODE)
  boost::optional<std::wstring>
  exercise_wstring(
    boost::optional<std::wstring> const& value)
  {
    return boost::optional<std::wstring>(*value + *value);
  }
#endif
}

BOOST_PYTHON_MODULE(boost_optional_ext)
{
  using boost_adaptbx::optional_conversions::to_and_from_python;
  to_and_from_python<boost::optional<bool> >();
  to_and_from_python<boost::optional<int> >();
  to_and_from_python<boost::optional<unsigned> >();
#if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
  to_and_from_python<boost::optional<std::size_t> >();
#endif
  to_and_from_python<boost::optional<float> >();
  to_and_from_python<boost::optional<double> >();
  to_and_from_python<boost::optional<std::string> >();
  boost::python::def("exercise", exercise);

#if !defined(BOOST_NO_STD_WSTRING) && defined(Py_USING_UNICODE)
  to_and_from_python<boost::optional<std::wstring> >();
  boost::python::def("exercise_wstring", exercise_wstring);
#endif
}
