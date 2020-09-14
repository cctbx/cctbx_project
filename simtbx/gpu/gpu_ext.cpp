#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <simtbx/gpu/structure_factors.h>

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

  } // namespace gpu

  BOOST_PYTHON_MODULE(simtbx_gpu_ext)
  {
    simtbx::gpu::structure_factor_wrapper::wrap();
  }
} // namespace simtbx
