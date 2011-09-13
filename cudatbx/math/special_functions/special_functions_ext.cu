#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cudatbx/math/special_functions/spherical_bessel_jn.h>

namespace cudatbx {
namespace math {
namespace special_functions {

namespace boost_python {

  void wrap_functions()
  {
    boost::python::def("cuda_spherical_bessel_jn",&cuda_spherical_bessel_jn);
  }

}}}

BOOST_PYTHON_MODULE(cudatbx_special_functions_ext)
{
  cudatbx::math::special_functions::boost_python::wrap_functions();
}
}
