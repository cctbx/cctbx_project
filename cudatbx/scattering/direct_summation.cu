#include <cudatbx/scattering/direct_summation.cuh>

namespace cudatbx {
namespace scattering {

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

    // reorganize coordinates
    int n_xyz = xyz.size();
    int padded_n_xyz = int(std::floor(n_xyz/padding + 1.0)) * padding;
    int size_xyz = 3 * padded_n_xyz;
    fType* h_xyz = new fType[size_xyz];
    for (int i=0; i<n_xyz; i++) {
      for (int j=0; j<3; j++) {
        h_xyz[j*padded_n_xyz + i] = fType(xyz[i][j]);
      }
    }

    // copy boundary layer solvent weights
    fType* h_solvent = new fType[padded_n_xyz];
    for (int i=0; i<n_xyz; i++) {
      h_solvent[i] = fType(solvent_weights[i]);
    }

    // reorganize h
    int n_h = h.size();
    int padded_n_h = int(std::floor(n_h/padding + 1.0)) * padding;
    int size_h = 3 * padded_n_h;
    fType* h_h = new fType[size_h];
    for (int i=0; i<n_h; i++) {
      for (int j=0; j<3; j++) {
        h_h[j*padded_n_h + i] = fType(h[i][j]);
      }
    }

    // reorganize rotations and translations
    // each rotation/translation pair is combined and padded to take up
    // 64 bytes so that a coalesced read will read two pairs
    int n_rt = translations.size();
    int size_rt = padded_size * n_rt;
    fType* h_rt = new fType[size_rt];
    for (int i=0; i<n_rt; i++) {
      for (int j=0; j<9; j++) {
        h_rt[padded_size*i + j] = fType(rotations[9*i + j]);
      }
      for (int j=0; j<3; j++) {
        h_rt[padded_size*i + j + 9] = fType(translations[i][j]);
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
    fType* h_a = new fType[f_size];
    fType* h_b = new fType[f_size];
    fType* h_c = new fType[n_types];
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
    fType* d_xyz;
    cudaSafeCall( cudaMalloc((void**)&d_xyz,size_xyz*sizeof(fType)) );
    cudaSafeCall( cudaMemcpy(d_xyz, h_xyz, size_xyz*sizeof(fType),
                             cudaMemcpyHostToDevice) );
    fType* d_solvent;
    cudaSafeCall( cudaMalloc((void**)&d_solvent,padded_n_xyz*sizeof(fType)) );
    cudaSafeCall( cudaMemcpy(d_solvent,h_solvent,padded_n_xyz*sizeof(fType),
                             cudaMemcpyHostToDevice) );
    fType* d_rt;
    cudaSafeCall( cudaMalloc((void**)&d_rt,size_rt*sizeof(fType)) );
    cudaSafeCall( cudaMemcpy(d_rt, h_rt, size_rt*sizeof(fType),
                             cudaMemcpyHostToDevice) );
    fType* d_h;
    cudaSafeCall( cudaMalloc((void**)&d_h,size_h*sizeof(fType)) );
    cudaSafeCall( cudaMemcpy(d_h, h_h, size_h*sizeof(fType),
                             cudaMemcpyHostToDevice) );
    // transfer data to constant memory
    // should combine d_n_types and d_n_terms into one transfer
    cudaSafeCall( cudaMemcpyToSymbol(d_a, h_a, f_size*sizeof(fType)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_b, h_b, f_size*sizeof(fType)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_c, h_c, n_types*sizeof(fType)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_n_types, &n_types, sizeof(int)) );
    cudaSafeCall( cudaMemcpyToSymbol(d_n_terms, &n_terms, sizeof(int)) );

    // allocate arrays for results if necessary
    if (sf_size == 0) {
      sf_size = n_h;
      cudaSafeCall( cudaMalloc((void**)&sf_real,n_h*sizeof(fType)) );
      cudaSafeCall( cudaMalloc((void**)&sf_imag,n_h*sizeof(fType)) );
    }
    else {
      assert(sf_size == n_h);
    }

    // run calculation
    int blocks_per_grid = (n_h + threads_per_block - 1)/threads_per_block;
    structure_factor_kernel<fType><<<blocks_per_grid,threads_per_block>>>
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
      fType* h_real = new fType[sf_size];
      fType* h_imag = new fType[sf_size];
      cudaSafeCall( cudaMemcpy(h_real,sf_real,sf_size*sizeof(fType),
                               cudaMemcpyDeviceToHost) );
      cudaSafeCall( cudaMemcpy(h_imag,sf_imag,sf_size*sizeof(fType),
                               cudaMemcpyDeviceToHost) );
      for (int i=0; i<sf_size; i++) {
        sf[i] = std::complex<double>(double(h_real[i]),double(h_imag[i]));
      }
      delete[] h_real;
      delete[] h_imag;
    }
    return sf;
  }

  /* ==========================================================================
   */

}
}
