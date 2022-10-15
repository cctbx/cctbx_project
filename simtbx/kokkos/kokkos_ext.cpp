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
      class_<simtbx::Kokkos::kokkos_instance>("gpu_instance",init<>() )
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
      class_<simtbx::Kokkos::kokkos_detector>("gpu_detector",init<>() )
        .def(init<int const&, const simtbx::nanoBragg::nanoBragg&>(
            ( arg("deviceId")=-1,arg("nanoBragg")),
             "Single panel constructor with data taken from nanoBragg instance\n"
             "deviceId is completely optional and ignored for Kokkos, rather\n"
             "the device is set by the gpu_instance class."))
        .def(init<int const&, dxtbx::model::Detector const &, dxtbx::model::Beam const &>(
            ( arg("deviceId")=-1,arg("detector"),arg("beam")),
             "Multipanel constructor with data taken from dxtbx objects\n"
             "deviceId is completely optional and ignored for Kokkos, rather\n"
             "the device is set by the gpu_instance class."))
        .def("show_summary",&simtbx::Kokkos::kokkos_detector::show_summary)
        .def("each_image_allocate",
              &simtbx::Kokkos::kokkos_detector::each_image_allocate,
             "Allocate large pixel arrays")
        .def("scale_in_place", &simtbx::Kokkos::kokkos_detector::scale_in_place,
             "Multiply by a scale factor on the GPU")
        .def("write_raw_pixels",&simtbx::Kokkos::kokkos_detector::write_raw_pixels,
             "Update raw_pixels on host with array from GPU")
        .def("get_raw_pixels",&simtbx::Kokkos::kokkos_detector::get_raw_pixels,
             "return multipanel detector raw pixels as a flex array")
        .def("get_whitelist_raw_pixels",
             (af::shared<double> (simtbx::Kokkos::kokkos_detector::*)(af::shared<std::size_t>))
             &simtbx::Kokkos::kokkos_detector::get_whitelist_raw_pixels,
            "return only those raw pixels requested by the whitelist selection, as a 1D flex array")
        .def("each_image_free", &simtbx::Kokkos::kokkos_detector::each_image_free)
        ;
    }
  };

  struct structure_factor_wrapper {
    static void
    wrap() {
      using namespace boost::python;
      class_<simtbx::Kokkos::kokkos_energy_channels>("gpu_energy_channels",init<>() )
        .def(init< const int& >(( arg("deviceId")=-1),
             "deviceId is optional and ignored for Kokkos, rather\n"
             "the device is set by the gpu_instance class."))
        .def("get_deviceID", &simtbx::Kokkos::kokkos_energy_channels::get_deviceID)
        .def("get_nchannels", &simtbx::Kokkos::kokkos_energy_channels::get_nchannels)
        .def("structure_factors_to_GPU_direct",
             &simtbx::Kokkos::kokkos_energy_channels::structure_factors_to_GPU_direct,
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
      using namespace simtbx::Kokkos;
      class_<simtbx::Kokkos::exascale_api>("exascale_api",no_init )
        .def(init<const simtbx::nanoBragg::nanoBragg&>(
            ( arg("nanoBragg"))))
        .def("allocate",&simtbx::Kokkos::exascale_api::allocate,
             "Allocate and transfer input data on the GPU")
        .def("add_energy_channel_from_gpu_amplitudes",
             &simtbx::Kokkos::exascale_api::add_energy_channel_from_gpu_amplitudes,
             (arg_("channel_number"), arg_("gpu_amplitudes"), arg_("gpu_detector"), arg_("weight")=1.0),
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spot contributions to the detector's accumulator array")
        .def("add_energy_channel_mask_allpanel",
             static_cast<void (exascale_api::*)(int const&,kokkos_energy_channels&,kokkos_detector&, af::shared<bool>) >
             (&exascale_api::add_energy_channel_mask_allpanel),
             (arg_("channel_number"), arg_("gpu_amplitudes"), arg_("gpu_detector"), arg_("pixel_active_mask_bools")),
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spots on mask==True pixels\n"
             "The pixel_active_mask_bools is a large array with one bool per detector pixel")
        .def("add_energy_channel_mask_allpanel",
             static_cast<void (exascale_api::*)(int const&,kokkos_energy_channels&,kokkos_detector&, af::shared<std::size_t> const) >
             (&exascale_api::add_energy_channel_mask_allpanel),
             (arg_("channel_number"), arg_("gpu_amplitudes"), arg_("gpu_detector"), arg_("pixel_active_list_ints")),
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spots on mask==True pixels\n"
             "The pixel_active_list_ints is a small array with integer-offset addresses for each pixel-of-interest")
        .def("add_energy_multichannel_mask_allpanel",
             static_cast<void (exascale_api::*)(af::shared<int> const,kokkos_energy_channels&,kokkos_detector&, af::shared<std::size_t> const,
             af::shared<double> const) >
             (&exascale_api::add_energy_multichannel_mask_allpanel),
             (arg_("ichannels"), arg_("gpu_amplitudes"), arg_("gpu_detector"), arg_("pixel_active_list_ints"), arg_("weights")),
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spots on mask==True pixels\n"
             "The pixel_active_list_ints is a small array with integer-offset addresses for each pixel-of-interest"
             "ichannels: for each nanoBragg source, the value instructs the simulation which channel in gpu_structure_factors"
             "to use for structure factor lookup.  If -1, skip this source wavelength."
             )
        .def("add_background", &simtbx::Kokkos::exascale_api::add_background,
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
