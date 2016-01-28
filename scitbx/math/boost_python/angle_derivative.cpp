#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../angle_derivative.h"

using namespace boost::python;

namespace scitbx { namespace math { namespace boost_python {

  void wrap_angle_derivative()
  {
    def("angle_derivative_wrt_vectors", &angle_derivative_wrt_vectors, (
      arg("u"),
      arg("v")));
  }

}}} // namespace scitbx::math::boost_python
