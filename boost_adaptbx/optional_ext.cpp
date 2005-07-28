#include <boost_adaptbx/optional_conversions.h>
#include <boost/python/module.hpp>

BOOST_PYTHON_MODULE(boost_optional_ext)
{
  using boost_adaptbx::optional_conversions::to_python;
  to_python<int>();
  to_python<unsigned>();
  to_python<float>();
  to_python<double>();
}
