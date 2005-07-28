#include <boost_adaptbx/optional_conversions.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

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
}

BOOST_PYTHON_MODULE(boost_optional_ext)
{
  using boost_adaptbx::optional_conversions::to_and_from_python;
  to_and_from_python<int>();
  to_and_from_python<unsigned>();
  to_and_from_python<float>();
  to_and_from_python<double>();
  boost::python::def("exercise", exercise);
}
