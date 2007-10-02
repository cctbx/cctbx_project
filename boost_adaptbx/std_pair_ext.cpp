#include <boost_adaptbx/std_pair_fwd.h>
#include <boost_adaptbx/optional_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost_adaptbx/std_pair_conversion.h>

namespace {
  std::pair<int, double> exercise(int i) {
    return std::make_pair(i*2, i/2.);
  }
}

BOOST_PYTHON_MODULE(std_pair_ext)
{
  using boost_adaptbx::std_pair_conversions::to_tuple;
  using boost::optional;

  to_tuple<int, double>();

  // next 4 needed for flex.find_partial_sum_xxx
  to_tuple< optional<std::size_t>, optional<int>      >();
  to_tuple< optional<std::size_t>, optional<unsigned> >();
  to_tuple< optional<std::size_t>, optional<float>    >();
  to_tuple< optional<std::size_t>, optional<double>   >();

  boost::python::def("exercise", exercise);
}
