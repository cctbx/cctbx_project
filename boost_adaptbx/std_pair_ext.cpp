#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost_adaptbx/std_pair_conversion.h>
#include <boost/optional.hpp>

namespace {
  std::pair<int, double> exercise(int i) {
    return std::make_pair(i*2, i/2.);
  }
}

BOOST_PYTHON_MODULE(std_pair_ext)
{
  using boost_adaptbx::std_pair_conversions::to_python;
  using boost::optional;

  to_python<int, double>();

  // next 4 needed for flex.find_partial_sum_xxx
  to_python< optional<std::size_t>, optional<int>      >();
  to_python< optional<std::size_t>, optional<unsigned> >();
  to_python< optional<std::size_t>, optional<float>    >();
  to_python< optional<std::size_t>, optional<double>   >();

  boost::python::def("exercise", exercise);
}
