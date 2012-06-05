#ifndef DIRECT_SUMMATION_CUH
#define DIRECT_SUMMATION_CUH

#include <cudatbx/scattering/direct_summation.h>
#include <cudatbx/cuda_base.cuh>

namespace cudatbx {
namespace scattering {

  __global__ void structure_factor_kernel
    (const int*, const float*, const float*, const int, const int,
     const float*, const int,
     const float*, const int,
     float*, float*);

  /* ==========================================================================
   */

}
}
#endif // DIRECT_SUMMATION_CUH
