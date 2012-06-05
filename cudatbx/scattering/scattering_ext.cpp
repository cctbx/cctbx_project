#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cudatbx/scattering/direct_summation.h>

namespace cudatbx { namespace scattering {

  struct direct_summation_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<cudatbx::scattering::direct_summation>("direct_summation",init<>() )
        .def("add",&cudatbx::scattering::direct_summation::add)
        .def("get_sum",&cudatbx::scattering::direct_summation::get_sum)
        ;
    }
  };

  }

  BOOST_PYTHON_MODULE(cudatbx_scattering_ext)
  {
    cudatbx::scattering::direct_summation_wrapper::wrap();
  }
}
