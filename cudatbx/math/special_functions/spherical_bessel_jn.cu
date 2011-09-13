#include <cudatbx/math/special_functions/spherical_bessel_jn.h>

// large values may result in too few registers per thread
const int threads_per_block = 256;

namespace cudatbx {
namespace math {
namespace special_functions {

  /* ==========================================================================
     Implementation of the spherical Bessel function based on the explicit
     formula for integer order, n
                                  n/2
                                  ___
                         1        \       k  a_2k(n + 0.5)
       j_n(z) = sin(z - --- n pi) /   (-1)  --------------- +
                         2        ---         z^(2k + 1)
                                  k=0

                                (n-1)/2
                                  ___
                         1        \       k  a_2k+1 (n + 0.5)
                cos(z - --- n pi) /   (-1)  ------------------
                         2        ---         z^(2k + 2)
                                  k=0

     where
                           (n + k)!
       a_k(n + 0.5) = ----------------- , k = 0, 1, ... , n
                       2^k k! (n - k)!

                    = 0                 , k = n + 1, n + 2, ...

     In the limit that z -> 0, the limiting form,

                 n
       j_n(z) = z  / (2n + 1)!!

     is used because of the denominator in the summations (z^(some power)).

     The transition between the limiting form and the explicit formula was
     determined empirically up to order 50.  The largest errors occur around
     the transition point, but should be lower than 1e-4 for high orders
     (lower errors for low orders).

     Notation is taken from the NIST Digital Library of Mathematical Functions
     (http://dlmf.nist.gov/), sections 10.49 (Equations 10.49.1 and 10.49.2)
     and 10.52 (Equation 10.52.1).

     The output is formatted into chunks of z, that is, all j0(z) are first,
     followed by all j1(z), then j2(z), etc.
     --------------------------------------------------------------------------
  */

  __global__ void spherical_bessel_jn_kernel
  (const int order, const double* z, const int n_z, double* j) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n_z) {
      double z_i = z[i];
      for (int n=0; n<order+1; n++) {
        j[n*n_z + i] = spherical_bessel_jn<int,double>(n,z_i);
      }
    }
  }

  __device__ double2 spherical_bessel_j0_j1(double z) {
    double sin_z, cos_z;
    sincos(z,&sin_z,&cos_z);
    z = 1.0/z;
    sin_z = sin_z*z;
    return make_double2(sin_z,sin_z*z - cos_z*z);
  }

  scitbx::af::shared<double> cuda_spherical_bessel_jn
  (const int& order, const scitbx::af::const_ref<double>& z, const int& gpu_id) {

    // start GPU
    cudaSafeCall( cudaSetDevice(gpu_id) );

    // allocate and initialize arrays
    int n_order = order + 1;
    int j_size = n_order * z.size();
    double * h_z, * h_j;
    double * d_z, * d_j;
    h_z = (double*)&z[0];
    h_j = new double[j_size];
    cudaSafeCall( cudaMalloc( (void**)&d_z, z.size() * sizeof(double) ) );
    cudaSafeCall( cudaMalloc( (void**)&d_j, j_size * sizeof(double) ) );
    cudaSafeCall( cudaMemcpy( d_z, h_z, z.size() * sizeof(double),
                              cudaMemcpyHostToDevice ) );

    // run kernel
    int blocks_per_grid = (z.size() + threads_per_block - 1)/threads_per_block;
    spherical_bessel_jn_kernel<<<blocks_per_grid,threads_per_block>>>
      (order,d_z,z.size(),d_j);

    // copy result from GPU
    cudaSafeCall( cudaMemcpy( h_j, d_j, j_size * sizeof(double),
                              cudaMemcpyDeviceToHost ) );
    scitbx::af::shared<double> jn((double*)&h_j[0], (double*)&h_j[0] + j_size);

    // clean up
    delete h_j;
    cudaSafeCall( cudaFree( d_z ) );
    cudaSafeCall( cudaFree( d_j ) );

    return jn;
  }

}
}
}
