#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost_adaptbx/std_pair_instantiation.h>

namespace {
  std::pair<int, double> exercise(int i) {
    return std::make_pair(i*2, i/2.);
  }
}

BOOST_PYTHON_MODULE(std_pair_ext)
{
  boost_adaptbx::std_pair_conversions::instantiate();

  boost::python::def("exercise", exercise);
}
