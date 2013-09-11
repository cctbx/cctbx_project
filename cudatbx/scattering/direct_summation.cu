#include <cudatbx/scattering/direct_summation.cuh>

namespace cudatbx {
namespace scattering {

  /* ==========================================================================
   */
  cudatbx::scattering::direct_summation::direct_summation() {
    // set host and device pointers to NULL
    h_xyz = NULL;
    h_solvent = NULL;
    h_h = NULL;
    h_rt = NULL;
    h_weights = NULL;
    h_scattering_type = NULL;
    h_a = NULL;
    h_b = NULL;
    h_c = NULL;

    d_xyz = NULL;
    d_solvent = NULL;
    d_h = NULL;
    d_rt = NULL;
    d_weights = NULL;
    d_scattering_type = NULL;

    amplitudes_allocated = false;
    h_real = NULL;
    h_imag = NULL;
    d_real = NULL;
    d_imag = NULL;

    workspace_allocated = false;
    d_workspace = NULL;
  }

  cudatbx::scattering::direct_summation::~direct_summation() {
    clear_arrays();
    clear_amplitudes();
    clear_workspace();
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_xyz
  (const scitbx::af::const_ref<scitbx::vec3<double> >& xyz) {
    // allocate memory if necessary
    if (n_xyz != xyz.size()) {
      clear_xyz();
      n_xyz = xyz.size();
      padded_n_xyz = cudatbx::calculate_padded_size(n_xyz,padding);
      size_xyz = 3 * padded_n_xyz;
      h_xyz = new fType[size_xyz];
      cudaSafeCall( cudaMalloc((void**)&d_xyz,size_xyz*sizeof(fType)) );
    }

    // convert values
    for (int i=0; i<n_xyz; i++) {
      for (int j=0; j<3; j++) {
        h_xyz[j*padded_n_xyz + i] = fType(xyz[i][j]);
      }
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_xyz, h_xyz, size_xyz*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_xyz() {
    delete[] h_xyz;
    cudaSafeCall( cudaFree(d_xyz) );
    h_xyz = NULL;
    d_xyz = NULL;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_solvent_weights
  (const scitbx::af::const_ref<double>& solvent_weights) {
    // allocate memory if necessary
    SCITBX_ASSERT (n_xyz == solvent_weights.size());
    if (n_solvent != solvent_weights.size()) {
      clear_solvent_weights();
      n_solvent = solvent_weights.size();
      h_solvent = new fType[padded_n_xyz];
      cudaSafeCall( cudaMalloc((void**)&d_solvent,padded_n_xyz*sizeof(fType)) );
    }

    // convert values
    for (int i=0; i<n_xyz; i++) {
      h_solvent[i] = fType(solvent_weights[i]);
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_solvent, h_solvent,
                             padded_n_xyz*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_solvent_weights() {
    delete[] h_solvent;
    cudaSafeCall( cudaFree(d_solvent) );
    h_solvent = NULL;
    d_solvent = NULL;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_hkl
  (const scitbx::af::const_ref<scitbx::vec3<double> >& h) {
    // allocate memory if necessary
    if (n_h != h.size()) {
      clear_hkl();
      n_h = h.size();
      padded_n_h = cudatbx::calculate_padded_size(n_h,padding);
      size_h = 3 * padded_n_h;
      h_h = new fType[size_h];
      cudaSafeCall( cudaMalloc((void**)&d_h,size_h*sizeof(fType)) );
    }

    // convert values
    for (int i=0; i<n_h; i++) {
      for (int j=0; j<3; j++) {
        h_h[j*padded_n_h + i] = fType(h[i][j]);
      }
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_h, h_h, size_h*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_hkl() {
    delete[] h_h;
    cudaSafeCall( cudaFree(d_h) );
    h_h = NULL;
    d_h = NULL;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_q
  (const scitbx::af::const_ref<double>& q) {
    // q data, use h variables
    if (n_h != q.size()) {
      clear_q();
      n_h = q.size();
      padded_n_h = cudatbx::calculate_padded_size(n_h,padding);
      size_h = padded_n_h;
      h_h = new fType[size_h];
      cudaSafeCall( cudaMalloc((void**)&d_h,size_h*sizeof(fType)) );
    }

    // convert values
    for (int i=0; i<n_h; i++) {
      h_h[i] = fType(q[i]);
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_h, h_h, size_h*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_q() {
    clear_hkl();
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_lattice
  (const scitbx::af::const_ref<double>& lattice_weights,
   const scitbx::af::const_ref<double>& lattice) {
    // lattice points, use rotation/translation
    if (n_rt != lattice_weights.size()) {
      clear_lattice();
      n_rt = lattice_weights.size();
      size_rt = cudatbx::calculate_padded_size(n_rt,padding);
      h_weights = new fType[size_rt];
      h_rt = new fType[3*size_rt];
      cudaSafeCall( cudaMalloc((void**)&d_weights,size_rt*sizeof(fType)) );
      cudaSafeCall( cudaMalloc((void**)&d_rt,3*size_rt*sizeof(fType)) );
    }

    // convert values
    for (int i=0; i<n_rt; i++) {
      h_weights[i] = fType(lattice_weights[i]/n_rt);
      for (int j=0; j<3; j++) {
        h_rt[j*size_rt + i] = fType(lattice[j*n_rt + i]);
      }
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_weights, h_weights, size_rt*sizeof(fType),
                             cudaMemcpyHostToDevice) );
    cudaSafeCall( cudaMemcpy(d_rt, h_rt, 3*size_rt*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_weights() {
    delete[] h_weights;
    cudaSafeCall( cudaFree(d_weights) );
    h_weights = NULL;
    d_weights = NULL;
  }

  void cudatbx::scattering::direct_summation::clear_lattice() {
    clear_weights();
    clear_rotations_translations();
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_rotations_translations
  (const scitbx::af::const_ref<double>& rotations,
   const scitbx::af::const_ref<scitbx::vec3<double> >& translations) {
    // each rotation/translation pair is combined and padded to take up
    // 64 bytes so that a coalesced read will read two pairs
    if (n_rt != translations.size()) {
      clear_rotations_translations();
      n_rt = translations.size();
      size_rt = padded_size * n_rt;
      h_rt = new fType[size_rt];
      cudaSafeCall( cudaMalloc((void**)&d_rt,size_rt*sizeof(fType)) );
    }

    // convert values
    for (int i=0; i<n_rt; i++) {
      for (int j=0; j<9; j++) {
        h_rt[padded_size*i + j] = fType(rotations[9*i + j]);
      }
      for (int j=0; j<3; j++) {
        h_rt[padded_size*i + j + 9] = fType(translations[i][j]);
      }
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_rt, h_rt, size_rt*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_rotations_translations() {
    delete[] h_rt;
    cudaSafeCall( cudaFree(d_rt) );
    h_rt = NULL;
    d_rt = NULL;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_scattering_types
  (const scitbx::af::const_ref<std::string>& scatterers,
   const cctbx::xray::scattering_type_registry& registry) {
    // allocate memory if necessary
    SCITBX_ASSERT (n_xyz == scatterers.size());
    if (n_scatterers != scatterers.size()) {
      clear_scattering_types();
      n_scatterers = scatterers.size();
      h_scattering_type = new int[padded_n_xyz];
      cudaSafeCall( cudaMalloc((void**)&d_scattering_type,
                               padded_n_xyz*sizeof(int)) );
    }

    // convert values
    for (int i=0; i<n_xyz; i++) {
      h_scattering_type[i] = registry.unique_index(scatterers[i]);
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpy(d_scattering_type,h_scattering_type,
                             padded_n_xyz*sizeof(int),cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_scattering_types() {
    delete[] h_scattering_type;
    cudaSafeCall( cudaFree(d_scattering_type) );
    h_scattering_type = NULL;
    d_scattering_type = NULL;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::set_scattering_type_registry
  (const cctbx::xray::scattering_type_registry& registry,
   const bool& complex_form_factor) {
    // convert form factors
    // add ordinary oxygen form factor at end for boundary layer solvent
    scitbx::af::shared<boost::optional
                       <cctbx::eltbx::xray_scattering::gaussian> >
      unique_gaussians = registry.unique_gaussians;
    n_types = unique_gaussians.size() + 1;
    n_terms = unique_gaussians[0].get().n_terms();
    f_size = n_types * n_terms;
    delete[] h_a;
    delete[] h_b;
    delete[] h_c;
    h_a = new fType[f_size];
    h_b = new fType[f_size];
    h_c = new fType[n_types];
    for (int i=0; i<f_size; i++) {
      h_a[i] = fType(0.0);
      h_b[i] = fType(0.0);
    }
    for (int i=0; i<n_types-1; i++) {
      for (int j=0; j<n_terms; j++) {
        h_a[i*n_terms + j] = fType(unique_gaussians[i].get().array_of_a()[j]);
        h_b[i*n_terms + j] = fType(unique_gaussians[i].get().array_of_b()[j]);
      }
      if (unique_gaussians[i].get().use_c()) {
        h_c[i] = fType(unique_gaussians[i].get().c());
      }
      else {
        h_c[i] = fType(0.0);
      }
    }

    // add form factor for boundary layer solvent
    cctbx::eltbx::xray_scattering::gaussian hoh =
      cctbx::eltbx::xray_scattering::wk1995("O",true).fetch();
    for (int i=0; i<hoh.array_of_a().size(); i++){
      h_a[(n_types-1)*n_terms + i] = fType(hoh.array_of_a()[i]);
      h_b[(n_types-1)*n_terms + i] = fType(hoh.array_of_b()[i]);
    }
    if (hoh.use_c()) {
      h_c[n_types-1] = fType(hoh.c());
    }
    else {
      h_c[n_types-1] = fType(0.0);
    }

    // transfer to GPU
    cudaSafeCall( cudaMemcpyToSymbol(dc_a, h_a, f_size*sizeof(fType)) );
    cudaSafeCall( cudaMemcpyToSymbol(dc_b, h_b, f_size*sizeof(fType)) );
    cudaSafeCall( cudaMemcpyToSymbol(dc_c, h_c, n_types*sizeof(fType)) );
    cudaSafeCall( cudaMemcpyToSymbol(dc_n_types, &n_types, sizeof(int)) );
    cudaSafeCall( cudaMemcpyToSymbol(dc_n_terms, &n_terms, sizeof(int)) );
    cudaSafeCall( cudaMemcpyToSymbol(dc_complex_form_factor,
                                     &complex_form_factor, sizeof(bool)) );
  }

  void cudatbx::scattering::direct_summation::clear_scattering_type_registry() {
    delete[] h_a;
    delete[] h_b;
    delete[] h_c;
    h_a = NULL;
    h_b = NULL;
    h_c = NULL;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::clear_arrays() {
    // clear pointers and set all pointers to NULL
    clear_xyz();
    clear_solvent_weights();
    clear_hkl();
    clear_rotations_translations();
    clear_weights();
    clear_scattering_types();
    clear_scattering_type_registry();
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::allocate_amplitudes() {
    if (amplitudes_allocated) {
      clear_amplitudes();
    }
    h_real = new fType[n_h];
    h_imag = new fType[n_h];
    cudaSafeCall( cudaMalloc((void**)&d_real,n_h*sizeof(fType)) );
    cudaSafeCall( cudaMalloc((void**)&d_imag,n_h*sizeof(fType)) );
    amplitudes_allocated = true;
  }

  void cudatbx::scattering::direct_summation::reset_amplitudes() {
    fType zero = fType(0.0);
    for (int i=0; i<n_h; i++) {
      h_real[i] = zero;
      h_imag[i] = zero;
    }
    cudaSafeCall( cudaMemcpy(d_real,h_real,n_h*sizeof(fType),
                             cudaMemcpyHostToDevice) );
    cudaSafeCall( cudaMemcpy(d_imag,h_imag,n_h*sizeof(fType),
                             cudaMemcpyHostToDevice) );
  }

  void cudatbx::scattering::direct_summation::clear_amplitudes() {
    delete[] h_real;
    delete[] h_imag;
    cudaSafeCall( cudaFree(d_real) );
    cudaSafeCall( cudaFree(d_imag) );
    h_real = NULL;
    h_imag = NULL;
    d_real = NULL;
    d_imag = NULL;
    amplitudes_allocated = false;
  }

  // --------------------------------------------------------------------------
  void cudatbx::scattering::direct_summation::allocate_workspace
  (const int& length) {
    if (workspace_allocated) {
      clear_workspace();
    }
    cudaSafeCall( cudaMalloc((void**)&d_workspace,length*sizeof(fType)) );
    workspace_allocated = true;
  }

  void cudatbx::scattering::direct_summation::clear_workspace() {
    cudaSafeCall( cudaFree(d_workspace) );
    d_workspace = NULL;
    workspace_allocated = false;
  }

  /* --------------------------------------------------------------------------
     reorganizes data and calls cuda
     padded to multiple of 128 bytes, (32 * sizeof(float or int))
  */
  void cudatbx::scattering::direct_summation::run_kernel() {
    int blocks_per_grid = cudatbx::calculate_blocks_per_grid(n_h,threads_per_block);
    structure_factor_kernel<fType><<<blocks_per_grid,threads_per_block>>>
      (d_scattering_type, d_xyz, d_solvent, n_xyz, padded_n_xyz,
       d_h, n_h, padded_n_h,
       d_rt, n_rt,
       d_real, d_imag);
  }

  void cudatbx::scattering::direct_summation::add
  (const scitbx::af::const_ref<std::string>& scatterers,
   const scitbx::af::const_ref<scitbx::vec3<double> >& xyz,
   const scitbx::af::const_ref<double>& solvent_weights,
   const scitbx::af::const_ref<scitbx::vec3<double> >& h,
   const scitbx::af::const_ref<double>& rotations,
   const scitbx::af::const_ref<scitbx::vec3<double> >& translations,
   const cctbx::xray::scattering_type_registry& registry,
   const bool& complex_form_factor) {

    // reorganize input data, allocates arrays, transfer to GPU, order matters
    set_xyz(xyz);
    set_solvent_weights(solvent_weights);
    set_hkl(h);
    set_rotations_translations(rotations,translations);
    set_scattering_types(scatterers,registry);
    set_scattering_type_registry(registry,complex_form_factor);

    // allocate arrays for results if necessary
    if (!amplitudes_allocated) {
      allocate_amplitudes();
      reset_amplitudes();
    }

    // run calculation
    run_kernel();

    // deallocate arrays
    clear_arrays();
  }

  /* --------------------------------------------------------------------------
     reorganizes data and calls cuda
     padded to multiple of 128 bytes, (32 * sizeof(float or int))

     "Rapid and accurate calculation of small-angle scattering profiles using
      the golden ratio"
     Watson, MC, Curtis, JE. J. Appl. Cryst. (2013). 46, 1171-1177

     solvent variables are used for weights and code is not optimal
     possibly subclass or split everything into functions
  */
  void cudatbx::scattering::direct_summation::prepare_saxs
  (const scitbx::af::const_ref<std::string>& scatterers,
   const scitbx::af::const_ref<scitbx::vec3<double> >& xyz,
   const scitbx::af::const_ref<double>& solvent_weights,
   const scitbx::af::const_ref<double>& q,
   const scitbx::af::const_ref<double>& lattice_weights,
   const scitbx::af::const_ref<double>& lattice,
   const cctbx::xray::scattering_type_registry& registry,
   const bool& complex_form_factor) {

    // reorganize input data, allocates arrays, transfer to GPU, order matters
    set_xyz(xyz);
    set_solvent_weights(solvent_weights);
    set_q(q);
    set_lattice(lattice_weights,lattice);
    set_scattering_types(scatterers,registry);
    set_scattering_type_registry(registry,complex_form_factor);

    // allocate arrays for results if necessary
    if (!amplitudes_allocated) {
      allocate_amplitudes();
    }
  }

  void cudatbx::scattering::direct_summation::run_saxs_kernel() {
    // allocate working space if necessary
    if (!workspace_allocated) {
      workspace_size = int(std::floor(n_h*n_rt/padding + 1.0)) * padding;
      allocate_workspace(3*workspace_size);
    }

    int blocks_per_grid = cudatbx::calculate_blocks_per_grid
      (n_rt,threads_per_block);
    expand_q_lattice_kernel<fType><<<blocks_per_grid,threads_per_block>>>
      (d_h, n_h,
       d_rt, n_rt, size_rt,
       d_workspace, workspace_size);
    blocks_per_grid = cudatbx::calculate_blocks_per_grid
      (n_h*n_rt,threads_per_block);
    saxs_kernel<fType><<<blocks_per_grid,threads_per_block>>>
      (d_scattering_type, d_xyz, d_solvent, n_xyz, padded_n_xyz,
       n_h, n_rt,
       d_workspace, workspace_size);
  }

  void cudatbx::scattering::direct_summation::run_solvent_saxs_kernel() {
    // allocate working space if necessary
    if (!workspace_allocated) {
      workspace_size = int(std::floor(n_h*n_rt/padding + 1.0)) * padding;
      allocate_workspace(7*workspace_size);
    }

    int blocks_per_grid = cudatbx::calculate_blocks_per_grid
      (n_rt,threads_per_block);
    expand_q_lattice_kernel<fType><<<blocks_per_grid,threads_per_block>>>
      (d_h, n_h,
       d_rt, n_rt, size_rt,
       d_workspace, workspace_size);
    blocks_per_grid = cudatbx::calculate_blocks_per_grid
      (n_h*n_rt,threads_per_block);
    solvent_saxs_kernel<fType><<<blocks_per_grid,threads_per_block>>>
      (d_scattering_type, d_xyz, d_solvent, n_xyz, padded_n_xyz,
       n_h, n_rt,
       d_workspace, workspace_size);
  }

  void cudatbx::scattering::direct_summation::run_collect_solvent_saxs_kernel
  (const double& c1, const double& c2) {
    assert(workspace_allocated);
    int blocks_per_grid = cudatbx::calculate_blocks_per_grid
      (n_h*n_rt,threads_per_block);
    collect_solvent_saxs_kernel<fType><<<blocks_per_grid,threads_per_block>>>
      (n_h, n_rt,fType(c1),fType(c2),d_workspace, workspace_size);
  }

  void cudatbx::scattering::direct_summation::sum_over_lattice() {
    int blocks_per_grid = cudatbx::calculate_blocks_per_grid
      (n_rt,threads_per_block);
    for (int i=0; i<n_h; i++) {
      cudatbx::math::weighted_sum_kernel<fType>
        <<<blocks_per_grid,threads_per_block,threads_per_block*sizeof(fType)>>>
        (&d_workspace[i*n_rt],d_weights,n_rt,&d_real[i]);
    }
  }

  /* --------------------------------------------------------------------------
     return total sum
  */
  scitbx::af::shared<std::complex<double> >
  cudatbx::scattering::direct_summation::get_sum() {
    scitbx::af::shared<std::complex<double> > sf(n_h);
    assert(amplitudes_allocated);
    cudaSafeCall( cudaMemcpy(h_real,d_real,n_h*sizeof(fType),
                             cudaMemcpyDeviceToHost) );
    cudaSafeCall( cudaMemcpy(h_imag,d_imag,n_h*sizeof(fType),
                             cudaMemcpyDeviceToHost) );
    for (int i=0; i<n_h; i++) {
      sf[i] = std::complex<double>(double(h_real[i]),double(h_imag[i]));
    }

    return sf;
  }

  /* ==========================================================================
   */

}
}
