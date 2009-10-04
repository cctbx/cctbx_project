#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost_adaptbx/optional_conversions.h>
#include <boost/optional.hpp>

namespace {

  boost::optional<int>
  exercise(
    boost::optional<double> const& value)
  {
    if (value) {
      if (*value == 13) return boost::optional<int>();
      return boost::optional<int>(static_cast<int>((*value)*3));
    }
    return boost::optional<int>(42);
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
  to_and_from_python<boost::optional<int> >();
  to_and_from_python<boost::optional<std::size_t> >();
  if (boost::python::converter::registry::query(
        boost::python::type_id<boost::optional<unsigned> >()) == 0) {
    to_and_from_python<boost::optional<unsigned> >();
  }
  to_and_from_python<boost::optional<float> >();
  to_and_from_python<boost::optional<double> >();
  to_and_from_python<boost::optional<std::string> >();
  boost::python::def("exercise", exercise);

#if !defined(BOOST_NO_STD_WSTRING) && defined(Py_USING_UNICODE)
  to_and_from_python<boost::optional<std::wstring> >();
  boost::python::def("exercise_wstring", exercise_wstring);
#endif
}
