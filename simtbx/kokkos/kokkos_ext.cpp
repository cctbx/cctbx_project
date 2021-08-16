#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include "simtbx/kokkos/kokkos_instance.h"
#include "simtbx/kokkos/detector.h"
#include "simtbx/kokkos/structure_factors.h"
#include "simtbx/kokkos/simulation.h"

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

  struct detector_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::Kokkos::kokkos_detector>("kokkos_detector",init<>() )
        .def(init<const simtbx::nanoBragg::nanoBragg&>(
            ( arg("nanoBragg"))))
//             "Single panel constructor with data taken from nanoBragg instance")
        .def(init<dxtbx::model::Detector const &, dxtbx::model::Beam const &>(
            ( arg("detector"),arg("beam"))))
//             "Multipanel constructor with data taken from dxtbx objects")
        .def("show_summary",&simtbx::Kokkos::kokkos_detector::show_summary)
        .def("each_image_allocate_cuda",
              &simtbx::Kokkos::kokkos_detector::each_image_allocate_cuda,
             "Allocate large pixel arrays")
        .def("scale_in_place_cuda", &simtbx::Kokkos::kokkos_detector::scale_in_place_cuda,
             "Multiply by a scale factor on the GPU")
        .def("write_raw_pixels_cuda",&simtbx::Kokkos::kokkos_detector::write_raw_pixels_cuda,
             "Update raw_pixels on host with array from GPU")
        .def("get_raw_pixels_cuda",&simtbx::Kokkos::kokkos_detector::get_raw_pixels_cuda,
             "return multipanel detector raw pixels as a flex array")
       .def("get_whitelist_raw_pixels_cuda",
             (af::shared<double> (simtbx::Kokkos::kokkos_detector::*)(af::shared<std::size_t>))
             &simtbx::Kokkos::kokkos_detector::get_whitelist_raw_pixels_cuda,
            "return only those raw pixels requested by the whitelist selection, as a 1D flex array")
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
        .def("print_Fhkl", &simtbx::Kokkos::kokkos_energy_channels::print_Fhkl,
             (arg_("channel"), arg_("first_element"), arg_("last_element"))
            )
        ;
    }
  };

  struct simulation_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::Kokkos::exascale_api>("exascale_api",no_init )
        .def(init<const simtbx::nanoBragg::nanoBragg&>(
            ( arg("nanoBragg"))))
        .def("allocate_cuda",&simtbx::Kokkos::exascale_api::allocate_cuda,
             "Allocate and transfer input data on the GPU")
        .def("add_energy_channel_from_kokkos_amplitudes_cuda",
             &simtbx::Kokkos::exascale_api::add_energy_channel_from_kokkos_amplitudes_cuda,
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spot contributions to the detector's accumulator array")
        .def("add_energy_channel_mask_allpanel_cuda",
             &simtbx::Kokkos::exascale_api::add_energy_channel_mask_allpanel_cuda,
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spots on mask==True pixels")
        .def("add_background_cuda", &simtbx::Kokkos::exascale_api::add_background_cuda,
             "Add a background field directly on the GPU")
        .def("show",&simtbx::Kokkos::exascale_api::show)
        ;
    }
  };

  } // namespace Kokkos

  BOOST_PYTHON_MODULE(simtbx_kokkos_ext) {
    simtbx::Kokkos::kokkos_instance_wrapper::wrap();
    simtbx::Kokkos::structure_factor_wrapper::wrap();
    simtbx::Kokkos::detector_wrapper::wrap();
    simtbx::Kokkos::simulation_wrapper::wrap();
  }
} // namespace simtbx
