#include "scitbx/array_family/boost_python/flex_fwd.h"
//#include <cudatbx/cuda_base.cuh>
#include "simtbx/kokkos/simulation.h"
#include "simtbx/kokkos/simulation_kernels.h"
#include "kokkostbx/kokkos_utils.h"
#include "scitbx/array_family/flex_types.h"

#define THREADS_PER_BLOCK_X 128
#define THREADS_PER_BLOCK_Y 1
#define THREADS_PER_BLOCK_TOTAL (THREADS_PER_BLOCK_X * THREADS_PER_BLOCK_Y)

namespace simtbx {
namespace Kokkos {

  namespace af = scitbx::af;
  //refactor later into helper file
/*  static cudaError_t cudaMemcpyVectorDoubleToDevice(CUDAREAL *dst, const double *src, size_t vector_items) {
        CUDAREAL * temp = new CUDAREAL[vector_items];
        for (size_t i = 0; i < vector_items; i++) {
                temp[i] = src[i];
        }
        cudaError_t ret = cudaMemcpy(dst, temp, sizeof(*dst) * vector_items, cudaMemcpyHostToDevice);
        delete temp;
        return ret;
  }*/

/*  exascale_api::~exascale_api(){
    cudaSafeCall(cudaSetDevice(SIM.device_Id));

        cudaSafeCall(cudaFree(cu_beam_vector));
        cudaSafeCall(cudaFree(cu_spindle_vector));
        cudaSafeCall(cudaFree(cu_source_X));
        cudaSafeCall(cudaFree(cu_source_Y));
        cudaSafeCall(cudaFree(cu_source_Z));
        cudaSafeCall(cudaFree(cu_source_I));
        cudaSafeCall(cudaFree(cu_source_lambda));
        cudaSafeCall(cudaFree(cu_a0));
        cudaSafeCall(cudaFree(cu_b0));
        cudaSafeCall(cudaFree(cu_c0));
        cudaSafeCall(cudaFree(cu_mosaic_umats));
        cudaSafeCall(cudaFree(cu_polar_vector));
  }
*/

} // Kokkos
} // simtbx
