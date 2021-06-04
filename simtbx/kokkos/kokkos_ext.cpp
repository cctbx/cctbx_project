#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include "simtbx/kokkos/kokkos_instance.h"

namespace simtbx { namespace Kokkos {

  struct kokkos_instance_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::Kokkos::kokkos_instance>("kokkos_instance",init<>() )
        .def(init< const int& >(( arg("deviceId"))))
        .def("finalize_kokkos", &simtbx::Kokkos::kokkos_instance::finalize_kokkos)
        ;
    }
  };

  } // namespace Kokkos

  BOOST_PYTHON_MODULE(simtbx_kokkos_ext) {
    simtbx::Kokkos::kokkos_instance_wrapper::wrap();
  }
} // namespace simtbx
