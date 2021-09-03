//#include "cudatbx/cuda_base.cuh"
#include "simtbx/kokkos/structure_factors.h"
#include "simtbx/kokkos/kokkos_utils.h"
#include <cstdio>

//using Kokkos::create_mirror_view;
//using Kokkos::deep_copy;
using Kokkos::parallel_for;
using Kokkos::RangePolicy;

namespace simtbx { namespace Kokkos {

  void
  kokkos_energy_channels::structure_factors_to_KOKKOS_detail(af::shared<double> linear_amplitudes){
    double * raw_ptr = linear_amplitudes.begin();

    vector_cudareal_t device_Fhkl( "device_Fhkl", linear_amplitudes.size() );
    transfer_shared2kokkos(device_Fhkl, linear_amplitudes);
    d_channel_Fhkl.push_back(device_Fhkl);

    if (d_channel_Fhkl.size()==1) { //first time through
      m_FhklParams = { h_range * k_range * l_range,
                             h_min, h_max, h_range, k_min, k_max, k_range, l_min, l_max, l_range };
    }
  }

  void
  kokkos_energy_channels::print_Fhkl(int channel, int first_element, int last_element) {
    auto view = d_channel_Fhkl[channel];
    parallel_for("print_Fhkl", RangePolicy<>(first_element, last_element), KOKKOS_LAMBDA (int i) {
      printf(" Fhkl[ %d ] = '%lf'\n", i, view(i));
    });
  }


} // Kokkos
} // simtbx
