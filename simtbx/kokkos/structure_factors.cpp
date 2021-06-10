//#include "cudatbx/cuda_base.cuh"
#include "simtbx/kokkos/structure_factors.h"
#include <cstdio>

using Kokkos::create_mirror_view;
using Kokkos::deep_copy;
using Kokkos::parallel_for;


namespace simtbx {
namespace Kokkos {

//  kokkos_energy_channels::kokkos_energy_channels(int const& deviceId){
//    h_deviceID = deviceId;
//    cudaSetDevice(deviceId);
//  }

  void
  kokkos_energy_channels::structure_factors_to_KOKKOS_detail(af::shared<double> linear_amplitudes){
    double * raw_ptr = linear_amplitudes.begin();

    view_vector_type device_Fhkl( "device_Fhkl", linear_amplitudes.size() );
    view_vector_type::HostMirror host_Fhkl = create_mirror_view( device_Fhkl );

    for (int i=0; i<linear_amplitudes.size(); ++i) {
      host_Fhkl(i) = (CUDAREAL) raw_ptr[i];
    }

    deep_copy(device_Fhkl, host_Fhkl);
    d_channel_Fhkl.push_back(device_Fhkl);

    if (d_channel_Fhkl.size()==1) { //first time through
      host_FhklParams = { h_range * k_range * l_range,
                             h_min, h_max, h_range, k_min, k_max, k_range, l_min, l_max, l_range };
    }
  }

  void
  kokkos_energy_channels::print_Fhkl(int channel, int elements) {
    auto view = d_channel_Fhkl[channel];
    parallel_for("print_Fhkl", elements, KOKKOS_LAMBDA (int i) {
      printf(" Fhkl[ %d ] = '%lf'\n", i, view(i));
    });
  }


} // Kokkos
} // simtbx
