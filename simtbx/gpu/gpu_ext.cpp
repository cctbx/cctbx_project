#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <simtbx/gpu/structure_factors.h>
#include <simtbx/gpu/detector.h>
#include <simtbx/gpu/simulation.h>

namespace simtbx { namespace gpu {

  struct structure_factor_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::gpu::gpu_energy_channels>("gpu_energy_channels",init<>() )
        .def(init< const int& >(( arg("deviceId"))))
        .def("get_deviceID",
             &simtbx::gpu::gpu_energy_channels::get_deviceID
            )
        .def("get_nchannels",
             &simtbx::gpu::gpu_energy_channels::get_nchannels
            )
        .def("structure_factors_to_GPU_direct_cuda",
             &simtbx::gpu::gpu_energy_channels::structure_factors_to_GPU_direct_cuda,
             (arg_("dummy_int"), arg_("indices"), arg_("amplitudes"))
            )
        ;
    }
  };

  struct detector_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      class_<simtbx::gpu::gpu_detector>("gpu_detector",init<>() )
        .def(init<int const&, const simtbx::nanoBragg::nanoBragg&>(
            ( arg("deviceId"),arg("nanoBragg"))))
//             "Single panel constructor with data taken from nanoBragg instance")
        .def(init<int const&, dxtbx::model::Detector const &, dxtbx::model::Beam const &>(
            ( arg("deviceId"),arg("detector"),arg("beam"))))
//             "Multipanel constructor with data taken from dxtbx objects")
        .def("get_deviceID", &simtbx::gpu::gpu_detector::get_deviceID
            )
        .def("show_summary",&simtbx::gpu::gpu_detector::show_summary)
        .def("each_image_allocate",
              &simtbx::gpu::gpu_detector::each_image_allocate,
             "Allocate large pixel arrays")
        .def("scale_in_place", &simtbx::gpu::gpu_detector::scale_in_place,
             "Apply a scale factor directly on the GPU")
        .def("write_raw_pixels",&simtbx::gpu::gpu_detector::write_raw_pixels,
             "Update raw_pixels on host with array from GPU")
        .def("get_raw_pixels",&simtbx::gpu::gpu_detector::get_raw_pixels,
             "return multipanel detector raw pixels as a flex array")
        .def("get_whitelist_raw_pixels",
              (af::shared<double> (simtbx::gpu::gpu_detector::*)(af::shared<std::size_t>))
              &simtbx::gpu::gpu_detector::get_whitelist_raw_pixels,
             "return only those raw pixels requested by the whitelist selection, as a 1D flex array")
        .def("each_image_free", &simtbx::gpu::gpu_detector::each_image_free)
        ;
    }
  };

  struct simulation_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      using namespace simtbx::gpu;
      class_<simtbx::gpu::exascale_api>("exascale_api",no_init )
        .def(init<const simtbx::nanoBragg::nanoBragg&>(
            ( arg("nanoBragg"))))
        .def("allocate_cuda",&simtbx::gpu::exascale_api::allocate_cuda,
             "Allocate and transfer input data on the GPU")
        .def("add_energy_channel_from_gpu_amplitudes_cuda",
             &simtbx::gpu::exascale_api::add_energy_channel_from_gpu_amplitudes_cuda,
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spot contributions to the detector's accumulator array")
        .def("add_energy_channel_mask_allpanel_cuda",
             static_cast<void (exascale_api::*)(int const&, gpu_energy_channels&, gpu_detector&, af::shared<int> const) >
             (&exascale_api::add_energy_channel_mask_allpanel_cuda),
             (arg_("channel_number"), arg_("gpu_amplitudes"), arg_("gpu_detector"), arg_("pixel_active_list_ints")),
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spots on mask==True pixels")
        .def("add_energy_channel_mask_allpanel_cuda",
             static_cast<void (exascale_api::*)(int const&, gpu_energy_channels&, gpu_detector&, af::shared<bool>) >
             (&exascale_api::add_energy_channel_mask_allpanel_cuda),
             (arg_("channel_number"), arg_("gpu_amplitudes"), arg_("gpu_detector"), arg_("pixel_active_mask_bools")),
             "Point to Fhkl at a new energy channel on the GPU, and accumulate Bragg spots on mask==True pixels")
        .def("add_background_cuda", &simtbx::gpu::exascale_api::add_background_cuda,
             (arg_("detector"), arg_("override_source")=-1),
             "Add a background field directly on the GPU")
        .def("show",&simtbx::gpu::exascale_api::show)
        ;
    }
  };

  } // namespace gpu

  BOOST_PYTHON_MODULE(simtbx_gpu_ext)
  {
    simtbx::gpu::structure_factor_wrapper::wrap();
    simtbx::gpu::detector_wrapper::wrap();
    simtbx::gpu::simulation_wrapper::wrap();
  }
} // namespace simtbx
