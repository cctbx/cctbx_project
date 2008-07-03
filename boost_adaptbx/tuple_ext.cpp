#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost_adaptbx/tuple_conversion.h>

namespace {

  boost::tuple<int, double> exercise(int i) {
    return boost::make_tuple(2*i, 0.5*i);
  }

}

BOOST_PYTHON_MODULE(boost_tuple_ext)
{
  using namespace boost_adaptbx::tuple_conversion;
  using namespace boost::python;
  to_python<boost::tuple<int, double> >();
  def("exercise", exercise);
}
