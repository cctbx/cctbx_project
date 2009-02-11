#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost_adaptbx/std_pair_conversion.h>
#include <boost_adaptbx/std_pair_fwd.h>

namespace {
  std::pair<int, double> exercise(std::pair<int, double> const &p) {
    return std::make_pair(p.first*2, p.second/2.);
  }
}

BOOST_PYTHON_MODULE(std_pair_ext)
{
  boost_adaptbx::std_pair_conversions::to_and_from_tuple<int, double>();
  boost::python::def("exercise", exercise);
}
