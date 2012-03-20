#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cudatbx/cuda_utilities.h>

namespace cudatbx {
namespace {

namespace boost_python {

  void wrap_functions()
  {
    boost::python::def("number_of_gpus",&number_of_gpus);
    boost::python::def("reset_gpu",&reset_gpu);
  }

}}

BOOST_PYTHON_MODULE(cudatbx_ext)
{
  cudatbx::boost_python::wrap_functions();
}
}
