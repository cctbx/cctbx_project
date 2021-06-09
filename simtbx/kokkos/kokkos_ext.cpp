#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include "simtbx/kokkos/kokkos_instance.h"
#include "simtbx/kokkos/structure_factors.h"

namespace simtbx { namespace Kokkos {

  struct kokkos_instance_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::Kokkos::kokkos_instance>("kokkos_instance",init<>() )
        .def(init< const int& >(( arg("deviceId"))))
        .def("get_deviceID", &simtbx::Kokkos::kokkos_instance::get_deviceID)
        .def("finalize_kokkos", &simtbx::Kokkos::kokkos_instance::finalize_kokkos)
        ;
    }
  };

  struct structure_factor_wrapper {
    static void
    wrap() {
      using namespace boost::python;
      class_<simtbx::Kokkos::kokkos_energy_channels>("kokkos_energy_channels",init<>() )
        .def("get_nchannels", &simtbx::Kokkos::kokkos_energy_channels::get_nchannels)
        .def("structure_factors_to_KOKKOS_direct_cuda",
             &simtbx::Kokkos::kokkos_energy_channels::structure_factors_to_KOKKOS_direct_cuda,
             (arg_("dummy_int"), arg_("indices"), arg_("amplitudes"))
            )
        ;
    }
  };

  } // namespace Kokkos

  BOOST_PYTHON_MODULE(simtbx_kokkos_ext) {
    simtbx::Kokkos::kokkos_instance_wrapper::wrap();
    simtbx::Kokkos::structure_factor_wrapper::wrap();
  }
} // namespace simtbx
