#include <cudatbx/scattering/direct_summation.cuh>

namespace cudatbx {
namespace scattering {

  const int threads_per_block = 1024;

  // form factors
  const int max_types = 50;
  const int max_terms = 10;
  __device__ __constant__ float d_a[max_types * max_terms];
  __device__ __constant__ float d_b[max_types * max_terms];
  __device__ __constant__ float d_c[max_types];
  __device__ __constant__ int d_n_types;
  __device__ __constant__ int d_n_terms;

  // constants
  __device__ __constant__ float two_pi = float(2.0)*CUDART_PI_F;
  const int padded_size = 16;
  __device__ __constant__ int d_padded_size = padded_size;

  /* ==========================================================================

     Memory properties for C2070, compute capability 2.0
     (CUDA C Programming Guide v 4.0, Appendix F)
     ---------------------------------------------------
     registers - register, r/w, "fast", per-multiprocessor   32K 4B registers
     local memory - r/w, "slow", per-thread                  512 KB
     __shared__ - shared memory, r/w, "fast", block-wide     48 KB
     __device__ - global memory, r/w, "slow", grid-wide,     6 GB
     __constant__ - constant memory, r, "fast", grid-wide    64 KB

     Shared memory is broken up into 32 banks and is interleaved into 32-bit
     words (4 bytes).  For example, an array of length 64 containing single
     precision values will have elements 0 and 32 in the same bank, 1 and 33
     in the same bank, etc.  To access shared memory with no conflicts (all
     threads get data with one read), each thread should read from a different
     bank, or have multiple threads read the same value in the same bank.  In
     the previous example, having all threads access element 0 or having each
     thread read a different element between 0 and 31, inclusive, will only
     require one read.  Accessing elements 0 and 32 will require two reads.

     Appendix F.4 describes the memory properties for compute capability 2.0
     devices in more detail and has figures for efficient memory access
     patterns.

     Basic approach
     --------------
     Each thread calculates the sum for one h, so each thread will
     independently loop over all atoms and put the sum into global memory

     All coordinates are loaded into global memory and then each thread copies
     sections into shared memory.  The kernel loops over all sections to sum
     over all atoms.  Rotation matrix/translation vectors pairs are also loaded
     and looped in the same manner.  Form factors are stored in constant memory

     The thread index is checked against the length of the array multiple times
     because all threads are used for reading atom data from global, but only
     threads whose index is less than the array length are needed for
     summation.  Additionaly, two __syncthreads() calls are required.  The
     first is to make sure all the atom data is copied into shared memory
     before any summation is started, and the second is to make sure all the
     summation is finished before the atom data in shared memory is replaced
     with new data.

     Data format for kernel
     ----------------------
     xyz = x_0 ... x_n y_0 ... y_n z_0 ... z_n
     solvent_weights = s_0 s_1 ... s_n
     h = h_0 ... h_n k_0 ... k_n l_0 ... l_n
     rt = r_00 ... r_08 t_00 ... t_02 ... r_10 ... t_n2
     a = a_00 a_01 a_02 ... a_n3 a_n4 a_n5
     b = ""
     c = c_0 c_1 ... c_n

     To facilitate coalesced reads from global memory, the data is grouped
     into sections.  For example, for xyz, all the x's come first, then all
     the y's, and lastly, all the z's.  When read from global memory, three
     coalesced reads will read in all the xyz's for a set of 32 atoms, one
     read from each section.  The size of the shared arrays is equal to the
     number of threads so that all threads will attempt to read from global
     memory.  There are checks against the actual length of available data.

     For the structue_factor_kernel, the general format of the loops is,

     -----------------------------
     | x_0 | x_1 | x_2 | x_3 | ...          xyz array in global memory
     -----------------------------
        |     |     |     |                 each thread stores one value into
        |     |     |     |                 shared memory
        V     V     V     V                 x[threadIdx.x] = xyz[current_atom];
     -----------------------------
     | x_0 | x_1 | x_2 | x_3 | ...          x array in shared memory
     -----------------------------
        |
        |-----|-----|-----|                 each thread reads one value
        V     V     V     V                 x_a = x[a];
     --------------------------------------------------------
     |each thread calculates its own sum with its registers |
     --------------------------------------------------------
        |     |     |     |
        |     |     |     |                 loop over all atoms
        V     V     V     V
     -----------------------------
     | r_0 | r_1 | r_2 | r_3 | ...          each thread copies its sums into
     -----------------------------          the structure factor arrays in
     -----------------------------          global memory
     | i_0 | i_1 | i_2 | i_3 | ...
     -----------------------------

     --------------------------------------------------------------------------
  */

  // kernel
  __global__ void structure_factor_kernel
  (const int* scattering_type, const float* xyz,
   const float* solvent_weights, const int n_xyz, const int padded_n_xyz,
   const float* h, const int n_h, const int padded_n_h,
   const float* rt, const int n_rt,
   float* sf_real, float* sf_imag) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    float h_i, k_i, l_i, stol_sq;
    float f[max_types];
    if (i < n_h) {
      // read h from global memory (stored in registers)
      h_i = h[               i];
      k_i = h[padded_n_h   + i];
      l_i = h[2*padded_n_h + i];

      // calculate form factors (stored in local memory)
      // last form factor is always for boundary solvent layer
      stol_sq = float(0.25) * (h_i*h_i + k_i*k_i + l_i*l_i);
      for (int type=0; type<d_n_types; type++) {
        f[type] = 0.0;
        for (int term=0; term<d_n_terms; term++) {
          f[type] += d_a[type*d_n_terms + term] *
            __expf(-d_b[type*d_n_terms + term] * stol_sq);
        }
        f[type] += d_c[type];
      }
    }

    // copy atoms into shared memory one chunk at a time and sum
    // all threads are used for reading data
    // shared arrays can be allocated at kernel invocation, but it requires
    // partitioning a big array (implement later)
    __shared__ float x[threads_per_block];
    __shared__ float y[threads_per_block];
    __shared__ float z[threads_per_block];
    __shared__ float solvent[threads_per_block];
    __shared__ int s_type[threads_per_block];
    __shared__ float rot_trans[threads_per_block];
    float real_sum = 0.0;
    float imag_sum = 0.0;
    float s,c,ff,xx,yy,zz,x_a,y_a,z_a;
    int current_atom, current_rt, rt_offset;

    for (int atom=0; atom<n_xyz; atom += blockDim.x) {
      current_atom = atom + threadIdx.x;
      // coalesce reads using threads, but don't read past n_xyz
      // one read for each variable should fill chunk of 32 atoms
      // total length = # of threads/block
      if (current_atom < n_xyz) {
        x[threadIdx.x] = xyz[                 current_atom];
        y[threadIdx.x] = xyz[padded_n_xyz   + current_atom];
        z[threadIdx.x] = xyz[2*padded_n_xyz + current_atom];
        solvent[threadIdx.x] = solvent_weights[current_atom];
        s_type[threadIdx.x] = scattering_type[current_atom];
      }

      // loop over all rotation/translation operators
      // one coalesced read will copy (# of threads)/(padded_size) rot/trans
      // since the number of threads is a multiple of 32, it will also always
      // be evenly divisible by padded_size
      for (int rt_i=0; rt_i<n_rt; rt_i += blockDim.x/d_padded_size) {
        current_rt = rt_i*d_padded_size + threadIdx.x;
        if (current_rt < n_rt*d_padded_size) {
          rot_trans[threadIdx.x] = rt[current_rt];
        }

        // wait for all data to be copied into shared memory
        __syncthreads();

        // then sum over all the atoms that are now available to all threads
        if (i < n_h) {
          for (int r=0; r<blockDim.x/d_padded_size; r++) {
            current_rt = rt_i + r;  // overall counter for rot/trans pairs
            if (current_rt < n_rt) {
              for (int a=0; a<blockDim.x; a++) {
                current_atom = atom + a;  // overall counter for atom number
                if (current_atom < n_xyz) {
                  x_a = x[a];  // transfer from shared memory to registers
                  y_a = y[a];  // might not be necessary due to cache
                  z_a = z[a];
                  rt_offset = r*d_padded_size;
                  // apply rotation and translation by expanding Rx + t
                  xx = (x_a*rot_trans[rt_offset    ] +
                        y_a*rot_trans[rt_offset + 1] +
                        z_a*rot_trans[rt_offset + 2] +
                        rot_trans[rt_offset + 9]);
                  yy = (x_a*rot_trans[rt_offset + 3] +
                        y_a*rot_trans[rt_offset + 4] +
                        z_a*rot_trans[rt_offset + 5] +
                        rot_trans[r*padded_size + 10]);;
                  zz = (x_a*rot_trans[rt_offset + 6] +
                        y_a*rot_trans[rt_offset + 7] +
                        z_a*rot_trans[rt_offset + 8] +
                        rot_trans[rt_offset + 11]);;
                  __sincosf(two_pi*(xx * h_i + yy * k_i + zz * l_i),&s,&c);
                  // bulk solvent correction in f
                  // boundary layer solvent scale in solvent
                  ff = f[s_type[a]] + solvent[a]*f[d_n_types-1];
                  real_sum += ff * c;
                  imag_sum += ff * s;
                }
              }
            }
          }
        }

        // wait before starting next chunk so data isn't changed for lagging threads
        __syncthreads();
      }
    }

    // transfer result to global memory
    if (i < n_h) {
      sf_real[i] = real_sum;
      sf_imag[i] = imag_sum;
    }
  }

  /* ==========================================================================
   */
  cudatbx::scattering::direct_summation::direct_summation() {
    sf_size = 0;
  }

  cudatbx::scattering::direct_summation::~direct_summation() {
    cudaSafeCall( cudaFree(sf_real) );
    cudaSafeCall( cudaFree(sf_imag) );
  }

  /* --------------------------------------------------------------------------
     reorganizes data and calls cuda
     padded to multiple of 128 bytes, (32 * sizeof(float or int))
  */
  void cudatbx::scattering::direct_summation::add
  (const scitbx::af::const_ref<std::string>& scatterers,
   const scitbx::af::const_ref<scitbx::vec3<double> >& xyz,
   const scitbx::af::const_ref<double>& solvent_weights,
   const scitbx::af::const_ref<scitbx::vec3<double> >& h,
   const scitbx::af::const_ref<double>& rotations,
   const scitbx::af::const_ref<scitbx::vec3<double> >& translations,
   const cctbx::xray::scattering_type_registry& registry) {

    int padding = 128 / sizeof(float);

    // reorganize coordinates
    int n_xyz = xyz.size();
    int padded_n_xyz = int(std::floor(n_xyz/padding + 1.0)) * padding;
    int size_xyz = 3 * padded_n_xyz;
    float* h_xyz = new float[size_xyz];
    for (int i=0; i<n_xyz; i++) {
      for (int j=0; j<3; j++) {
        h_xyz[j*padded_n_xyz + i] = float(xyz[i][j]);
      }
    }

    // copy boundary layer solvent weights
    float* h_solvent = new float[padded_n_xyz];
    for (int i=0; i<n_xyz; i++) {
      h_solvent[i] = float(solvent_weights[i]);
    }

    // reorganize h
    int n_h = h.size();
    int padded_n_h = int(std::floor(n_h/padding + 1.0)) * padding;
    int size_h = 3 * padded_n_h;
    float* h_h = new float[size_h];
    for (int i=0; i<n_h; i++) {
      for (int j=0; j<3; j++) {
        h_h[j*padded_n_h + i] = float(h[i][j]);
      }
    }

    // reorganize rotations and translations
    // each rotation/translation pair is combined and padded to take up
    // 64 bytes so that a coalesced read will read two pairs
    int n_rt = translations.size();
    int size_rt = padded_size * n_rt;
    float* h_rt = new float[size_rt];
    for (int i=0; i<n_rt; i++) {
      for (int j=0; j<9; j++) {
        h_rt[padded_size*i + j] = float(rotations[9*i + j]);
      }
      for (int j=0; j<3; j++) {
        h_rt[padded_size*i + j + 9] = float(translations[i][j]);
      }
    }

    // convert scattering types and form factors
    // add ordinary oxygen form factor at end for boundary layer solvent
    int* h_scattering_type = new int[padded_n_xyz];
    for (int i=0; i<n_xyz; i++) {
      h_scattering_type[i] = registry.unique_index(scatterers[i]);
    }
    scitbx::af::shared<boost::optional
                       <cctbx::eltbx::xray_scattering::gaussian> >
      unique_gaussians = registry.unique_gaussians;
    int n_types = unique_gaussians.size() + 1;
    int n_terms = unique_gaussians[0].get().n_terms();
    int f_size = n_types * n_terms;
    float* h_a = new float[f_size];
    float* h_b = new float[f_size];
    float* h_c = new float[n_types];
    for (int i=0; i<f_size; i++) {
      h_a[i] = 0.0;
      h_b[i] = 0.0;
    }
    for (int i=0; i<n_types-1; i++) {
      for (int j=0; j<n_terms; j++) {
        h_a[i*n_terms + j] = unique_gaussians[i].get().array_of_a()[j];
        h_b[i*n_terms + j] = unique_gaussians[i].get().array_of_b()[j];
      }
      if (unique_gaussians[i].get().use_c()) {
        h_c[i] = unique_gaussians[i].get().c();
      }
      else {
        h_c[i] = float(0.0);
      }
    }
    // add form factor for boundary layer solvent (# of terms may be different)
    cctbx::eltbx::xray_scattering::gaussian hoh =
      cctbx::eltbx::xray_scattering::wk1995("O",true).fetch();
    for (int i=0; i<hoh.array_of_a().size(); i++){
      h_a[(n_types-1)*n_terms + i] = hoh.array_of_a()[i];
      h_b[(n_types-1)*n_terms + i] = hoh.array_of_b()[i];
    }
    if (hoh.use_c()) {
      h_c[n_types-1] = hoh.c();
    }
    else {
      h_c[n_types-1] = float(0.0);
    }

    // transfer data to global memory
    int* d_scattering_type;
    cudaSafeCall( cudaMalloc((void**)&d_scattering_type,padded_n_xyz*sizeof(int)) );
    cudaSafeCall( cudaMemcpy(d_scattering_type,h_scattering_type,
                             padded_n_xyz*sizeof(int),cudaMemcpyHostToDevice) );
    float* d_xyz;
    cudaSafeCall( cudaMalloc((void**)&d_xyz,size_xyz*sizeof(float)) );
    cudaSafeCall( cudaMemcpy(d_xyz, h_xyz, size_xyz*sizeof(float),
                             cudaMemcpyHostToDevice) );
    float* d_solvent;
    cudaSafeCall( cudaMalloc((void**)&d_solvent,padded_n_xyz*sizeof(float)) );
    cudaSafeCall( cudaMemcpy(d_solvent,h_solvent,padded_n_xyz*sizeof(float),
                             cudaMemcpyHostToDevice) );
    float* d_rt;
    cudaSafeCall( cudaMalloc((void**)&d_rt,size_rt*sizeof(float)) );
    cudaSafeCall( cudaMemcpy(d_rt, h_rt, size_rt*sizeof(float),
                             cudaMemcpyHostToDevice) );
    float* d_h;
    cudaSafeCall( cudaMalloc((void**)&d_h,size_h*sizeof(float)) );
    cudaSafeCall( cudaMemcpy(d_h, h_h, size_h*sizeof(float),
                             cudaMemcpyHostToDevice) );
    // transfer data to constant memory
    // should combine d_n_types and d_n_terms into one transfer
    cudaSafeCall( cudaMemcpyToSymbol(d_a, h_a, f_size*sizeof(float)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_b, h_b, f_size*sizeof(float)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_c, h_c, n_types*sizeof(float)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_n_types, &n_types, sizeof(int)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_n_terms, &n_terms, sizeof(int)) );

    // allocate arrays for results if necessary
    if (sf_size == 0) {
      sf_size = n_h;
      cudaSafeCall( cudaMalloc((void**)&sf_real,n_h*sizeof(float)) );
      cudaSafeCall( cudaMalloc((void**)&sf_imag,n_h*sizeof(float)) );
    }
    else {
      assert(sf_size == n_h);
    }

    // run calculation
    int blocks_per_grid = (n_h + threads_per_block - 1)/threads_per_block;
    structure_factor_kernel<<<blocks_per_grid,threads_per_block>>>
      (d_scattering_type, d_xyz, d_solvent, n_xyz, padded_n_xyz,
       d_h, n_h, padded_n_h,
       d_rt, n_rt,
       sf_real, sf_imag);

    // clean up
    delete[] h_xyz;
    delete[] h_solvent;
    delete[] h_h;
    delete[] h_rt;
    delete[] h_scattering_type;
    delete[] h_a;
    delete[] h_b;
    delete[] h_c;
    cudaSafeCall( cudaFree(d_h) );
    cudaSafeCall( cudaFree(d_xyz) );
    cudaSafeCall( cudaFree(d_solvent) );
    cudaSafeCall( cudaFree(d_rt) );
    cudaSafeCall( cudaFree(d_scattering_type) );
  }

  /* --------------------------------------------------------------------------
     return total sum
  */
  scitbx::af::shared<std::complex<double> >
  cudatbx::scattering::direct_summation::get_sum() {
    scitbx::af::shared<std::complex<double> > sf(sf_size);
    if (sf_size != 0) {
      float* h_real = new float[sf_size];
      float* h_imag = new float[sf_size];
      cudaSafeCall( cudaMemcpy(h_real,sf_real,sf_size*sizeof(float),
                               cudaMemcpyDeviceToHost) );
      cudaSafeCall( cudaMemcpy(h_imag,sf_imag,sf_size*sizeof(float),
                               cudaMemcpyDeviceToHost) );
      for (int i=0; i<sf_size; i++) {
        sf[i] = std::complex<double>(double(h_real[i]),double(h_imag[i]));
      }
      delete[] h_real;
      delete[] h_imag;
    }
    return sf;
  }

}
}
