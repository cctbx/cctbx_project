import copy,math
from scitbx.array_family import flex
from cctbx.xray import ext
from libtbx.utils import Sorry

# alternative implementations of direct summation using parallel architecture
# example usage is given in cctbx/regression/tst_direct_scaling

class pprocess:
  """
  Software requirements for running direct summation using Cuda parallel processing units:
  1. numpy package (http://numpy.scipy.org) version 1.5.1 OK
  2. pycuda package, version 2011.1 OK
        more information:  http://wiki.tiker.net/PyCuda/Installation/Linux
                           http://documen.tician.de/pycuda/tutorial.html
        suggested install scheme for Linux (not necessarily exact):
          ../base/bin/python configure.py --cuda-root=/usr/common/usg/cuda/3.2 \
            --cudadrv-lib-dir=/usr/common/usg/nvidia-driver-util/3.2/lib64 \
            --boost-inc-dir=$HOME/boostbuild/include \
            --boost-lib-dir=$HOME/boostbuild/lib \
            --boost-python-libname=boost_python \
            --boost-thread-libname=boost_thread
          make install
  3.  gcc 4.4.2 or higher is required for Linux build of pycuda 2011.1
  4.  boost.python is required for pycuda; cctbx-installed version is probably OK, not tested.
        tests were performed with separately-installed boost:
          cd boost_1_45_0
          ./bootstrap.sh --prefix=$HOME/boostbuild --libdir=$HOME/boostbuild/lib --with-python=$HOME/build/base/bin/python --with-libraries=signals,thread,python
          ./bjam variant=release link=shared install
  5.  cuda 3.2 separately installed is required for pycuda 2011.1
        (http://developer.nvidia.com/object/cuda_3_2_downloads.html)

  Suggested hardware:
    as tested, the Nvidia Tesla C2050 (Fermi, compute capability 2.0)
               gave 50-fold performance improvement over CPU.
    The Nvidia Tesla C1060 (compute capability 1.3) gave 15-fold performance improvement.
  """

  def print_diagnostics(self):
    for scatterer in self.scatterers:
      print scatterer
      print scatterer.fp
      print scatterer.fdp
    self.registry.show()
    d_star_sq = 1.
    uniqueff = self.registry.unique_form_factors_at_d_star_sq(d_star_sq)
    uniqueix = self.registry.unique_indices(self.scatterers)
    for x in xrange(len(self.scatterers)):
      print self.scatterers[x]
      print self.scatterers[x].scattering_type
      print uniqueff[uniqueix[x]]
      constant_term = self.registry.gaussian(self.scatterers[x].scattering_type).c()
      print "constant",constant_term
      parameters = self.registry.gaussian(self.scatterers[x].scattering_type).parameters()
      x_sq = d_star_sq / 4.;
      unique_form_factor_at_x_sq = copy.copy(constant_term)
      for x in xrange(0,len(parameters),2):
        unique_form_factor_at_x_sq += parameters[x] * math.exp(-parameters[x+1]*x_sq)
      print unique_form_factor_at_x_sq

  def prepare_miller_arrays_for_cuda(self,numpy,verbose=False):

    """ The number of miller indices and atoms each must be an exact
        multiple of the BLOCKSIZE (32), so zero-padding is employed"""

    # Transfer data from flex arrays to numpy & zero-pad
    # Marshal the miller indices into numpy arrays for use in CUDA

    # Miller indices
    flat_mix = self.miller_indices.as_vec3_double().as_double().as_numpy_array()
    if self.miller_indices.size()%FHKL_BLOCKSIZE > 0:
      newsize = 3*FHKL_BLOCKSIZE* (1+(self.miller_indices.size()/FHKL_BLOCKSIZE))
      if verbose: print "Miller indices: resetting %s array from size %d to %d"%(
                        flat_mix.dtype,flat_mix.shape[0],newsize)
      flat_mix.resize(( newsize,))

    self.n_flat_hkl = flat_mix.shape[0]/3
    self.flat_mix = flat_mix

  def prepare_scattering_sites_for_cuda(self,numpy,verbose=False):

    # Transfer data from flex arrays to numpy & zero-pad
    # Marshal the fractional coordinates, and weights
    # into numpy arrays for use in CUDA

    # inspect the unique scatterers in the registry, determine number of electrons
    n_electrons = self.registry.unique_form_factors_at_d_star_sq(0.)

    # scatterer id, listed in increasing order of electron content
    self.scatterers.increasing_order = flex.sort_permutation(n_electrons)
    self.scatterers.number_of_types = len(self.scatterers.increasing_order)
    self.scatterers.unsorted_counts = []
    self.scatterers.sorted_ranges = []

    # Unique type labels for the scatterers
    uniqueix = self.registry.unique_indices(self.scatterers)
    for s in xrange(self.scatterers.number_of_types):
      self.scatterers.unsorted_counts.append( uniqueix.count(s) )
      self.scatterers.sorted_ranges.append([0,0])
    for s in xrange(self.scatterers.number_of_types):
      self.scatterers.sorted_ranges[self.scatterers.increasing_order[s]][1]=\
        self.scatterers.sorted_ranges[self.scatterers.increasing_order[s]][0]+\
        self.scatterers.unsorted_counts[s]
      if s+1 < self.scatterers.number_of_types:
        self.scatterers.sorted_ranges[self.scatterers.increasing_order[s+1]][0]=\
        self.scatterers.sorted_ranges[self.scatterers.increasing_order[s]][1]

    if verbose:
      print list(self.scatterers.increasing_order)
      print self.scatterers.number_of_types
      print self.scatterers.unsorted_counts
      print self.scatterers.sorted_ranges

    uniqueix_sort_order = flex.sort_permutation(uniqueix)
    sorted_uniqueix = uniqueix.select(uniqueix_sort_order
      ).as_numpy_array().astype(numpy.uint32)

    # Scattering sites
    sites = self.scatterers.extract_sites()
    n_sites = sites.size()
    sorted_flat_sites = sites.select(uniqueix_sort_order).as_double().as_numpy_array()
    #print sorted_flat_sites.dtype,sorted_flat_sites.shape

    # weights: weight = site_multiplicity * occupancy
    # for example, an atom on a three-fold has site_multiplicity 1/3
    sorted_weights = flex.double([S.weight() for S in self.scatterers]).select(
      uniqueix_sort_order).as_numpy_array()

    sorted_u_iso = self.scatterers.extract_u_iso().select(
      uniqueix_sort_order).as_numpy_array()

    if n_sites%FHKL_BLOCKSIZE > 0:
      newsize = FHKL_BLOCKSIZE* (1+(n_sites/FHKL_BLOCKSIZE))
      if verbose:
        print "Scatterer xyzs: resetting %s array from size %d to %d"%(
              sorted_flat_sites.dtype,sorted_flat_sites.shape[0],3*newsize)
        print "weights   : resetting %s array from size %d to %d"%(
              sorted_weights.dtype,sorted_weights.shape[0],newsize)
        print
      sorted_flat_sites.resize(( 3*newsize,))
      # zero-padding for weights guarantees no effect from padded sites
      sorted_weights.resize(( newsize,))
      sorted_u_iso.resize(( newsize,))
      sorted_uniqueix.resize((newsize,))

    self.n_sites = n_sites
    self.flat_sites = sorted_flat_sites
    self.weights = sorted_weights
    self.u_iso = sorted_u_iso
    self.uniqueix = sorted_uniqueix

  def prepare_gaussians_symmetries_cell(self,numpy):
    # Marshal the scattering types, form factors, symmetry elements,
    # and metrical matrix into numpy arrays for use in CUDA.

    # Gaussian expansion for the unique scattering types
    gaussians = flex.double()

    self.n_gaussians = 0
    self.n_terms_in_sum = 0
    for gaussian in self.registry.unique_gaussians_as_list():
      self.n_gaussians += 1
      #constant_term
      gaussians.append(gaussian.c())
      #non-constant terms in sum
      terms = gaussian.parameters()
      if self.n_terms_in_sum != 0:
        assert len(terms) == 2 * self.n_terms_in_sum
        assert len(terms) % 2 == 0
      self.n_terms_in_sum = len(terms)/2
      for item in terms:
        gaussians.append(item)

    self.gaussians = gaussians.as_numpy_array()
    self.g_stride = 1 + 2 *self.n_terms_in_sum
    assert self.gaussians.shape[0] == self.n_gaussians * self.g_stride

    # Space group symmetry
    self.nsymop = self.space_group.order_z()
    self.sym_stride = 12 # numbers to uniquely specify the symop
    symmetry = flex.double()

    for symop in self.space_group:
      for item in symop.r().as_double():  symmetry.append(item)
      for item in symop.t().as_double():  symmetry.append(item)
    self.symmetry = symmetry.as_numpy_array()
    assert self.symmetry.shape[0] == self.nsymop * self.sym_stride

    # Unit cell dimensions --> metrical matrix --> establishes dstar for given hkl
    self.reciprocal_metrical_matrix = numpy.array(
      self.unit_cell.reciprocal_metrical_matrix() )

  def validate_the_inputs(self,manager,cuda):

    # Assess limits due to the current implementation:
    # no support for ADP's, f' or f"

    PyNX_implementation = "Fhkl = Sum(atoms) (weight * exp(2*pi*i(H*X)))"
    # reference http://pynx.sourceforge.net (http://arxiv.org/abs/1010.2641v1)

    Present_implementation = """Fhkl = Sum(atoms(X)) (f0(d*) * weight *
        Sum(symops(SR,ST)) exp(2*pi*i(H*SR*X+H*ST) + DW-factor)
      )"""

    Future_implementation = """Fhkl = Sum(atoms(X)) ((f0(d*)+f'+i*f") * weight *
        Sum(symops(SR,ST)) exp(2*pi*i(H*SR*X+H*ST) + DW-factor)
      )"""
    try:
      assert len(self.miller_indices) > 0
    except:
      raise Sorry("There must be at least one Miller index""")

    self.use_debye_waller = (
      self.scatterers.extract_u_iso().count(0.0)<len(self.scatterers) )

    try:
      assert self.scatterers.extract_use_u_aniso().count(True)==0
    except:
      raise Sorry("As presently implemented parallel processor direct summation doesn't support anisotropic displacement factors""")

    try:
      assert [S.fp for S in self.scatterers].count(0.0)==len(self.scatterers)
      assert [S.fdp for S in self.scatterers].count(0.0)==len(self.scatterers)
    except:
      raise Sorry("As presently implemented parallel processor direct summation doesn't support anomalous scatterers""")

    # Assess limits based on global memory size of parallel unit
    n_atoms = len(self.scatterers)
    if n_atoms%FHKL_BLOCKSIZE > 0:
       n_atoms = FHKL_BLOCKSIZE* (1+(n_atoms/FHKL_BLOCKSIZE))
    n_hkl = len(self.miller_indices)
    if n_hkl%FHKL_BLOCKSIZE > 0:
       n_hkl = FHKL_BLOCKSIZE* (1+(n_hkl/FHKL_BLOCKSIZE))
    global_memory_atoms = 36 * n_atoms # refers to mod_fhkl_str CUDA code
    global_memory_hkl = 40 * n_hkl # refers to mod_fhkl_str CUDA code

    totalmem = cuda.Device(0).total_memory()
    if 1.2*(global_memory_atoms + global_memory_hkl) > totalmem:
      raise Sorry("Atoms use %d bytes; hkl indices use %d bytes, exceeding total GPU memory %d."%(
        global_memory_atoms,global_memory_hkl,totalmem))

    # Assess limit based on total available threads
    totalthreads = cuda.Device(0).get_attribute(
                     cuda.device_attribute.MAX_BLOCK_DIM_X) * \
                   cuda.Device(0).get_attribute(
                     cuda.device_attribute.MAX_GRID_DIM_X)
    if n_hkl > totalthreads:
      raise Sorry("""%d reflection values will exceed the maximum number of
        threads for a single kernel invocation, %d"""%(
         n_hkl,totalthreads))

  def __init__(self,instance,algorithm,verbose=False):
    """
    instance: an instance of class from_scatterers_direct(cctbx.xray.structure_factors.manager.managed_calculation_base)
    algorithm: an instance of class direct_summation_cuda_platform(direct_summation_simple)
    """

    self.scatterers = instance._xray_structure.scatterers()
    self.registry = instance._xray_structure.scattering_type_registry()
    self.miller_indices = instance._miller_set.indices()
    self.unit_cell = instance._miller_set.unit_cell()
    self.space_group = instance._miller_set.space_group()

    if verbose: self.print_diagnostics() # some diagnostics used for development

    if hasattr(algorithm,"simple"):
      instance._results = ext.structure_factors_simple(
      self.unit_cell,
      instance._miller_set.space_group(),
      self.miller_indices,
      self.scatterers,
      self.registry); return

    if hasattr(algorithm,"pycuda"):
      import pycuda.driver as cuda
      from pycuda.compiler import SourceModule

      self.validate_the_inputs(instance,cuda)

      self.prepare_miller_arrays_for_cuda(algorithm.numpy)

      self.prepare_scattering_sites_for_cuda(algorithm.numpy)

      self.prepare_gaussians_symmetries_cell(algorithm.numpy)

      assert cuda.Device.count() >= 1

      device = cuda.Device(0)
      WARPSIZE=device.get_attribute(cuda.device_attribute.WARP_SIZE) # 32
      MULTIPROCESSOR_COUNT=device.get_attribute(cuda.device_attribute.MULTIPROCESSOR_COUNT)

      sort_mod = SourceModule(mod_fhkl_sorted%(self.gaussians.shape[0],
                         self.symmetry.shape[0],
                         self.sym_stride,
                         self.g_stride,
                         int(self.use_debye_waller),
                         self.nsymop),
                         )

      r_m_m_address = sort_mod.get_global("reciprocal_metrical_matrix")[0]
      cuda.memcpy_htod(r_m_m_address, self.reciprocal_metrical_matrix)

      gaussian_address = sort_mod.get_global("gaussians")[0]
      cuda.memcpy_htod(gaussian_address, self.gaussians)

      symmetry_address = sort_mod.get_global("symmetry")[0]
      cuda.memcpy_htod(symmetry_address, self.symmetry)

      CUDA_fhkl = sort_mod.get_function("CUDA_fhkl")

      intermediate_real = algorithm.numpy.zeros((self.n_flat_hkl,),algorithm.numpy.float64)
      intermediate_imag = algorithm.numpy.zeros((self.n_flat_hkl,),algorithm.numpy.float64)
      for x in xrange(self.scatterers.number_of_types):
        fhkl_real = algorithm.numpy.zeros((self.n_flat_hkl,),algorithm.numpy.float64)
        fhkl_imag = algorithm.numpy.zeros((self.n_flat_hkl,),algorithm.numpy.float64)

        CUDA_fhkl(cuda.InOut(fhkl_real),
                 cuda.InOut(fhkl_imag),
                 cuda.In(self.flat_sites),
                 cuda.In(self.weights),
                 cuda.In(self.u_iso),
                 algorithm.numpy.uint32(self.scatterers.increasing_order[x]),
                 algorithm.numpy.uint32(self.scatterers.sorted_ranges[x][0]),
                 algorithm.numpy.uint32(self.scatterers.sorted_ranges[x][1]),
                 cuda.In(self.flat_mix),
                 block=(FHKL_BLOCKSIZE,1,1),
                 grid=((self.n_flat_hkl/FHKL_BLOCKSIZE,1)))

        intermediate_real += fhkl_real
        intermediate_imag += fhkl_imag

      flex_fhkl_real = flex.double(intermediate_real[0:len(self.miller_indices)])
      flex_fhkl_imag = flex.double(intermediate_imag[0:len(self.miller_indices)])

      instance._results = fcalc_container(flex.complex_double(flex_fhkl_real,flex_fhkl_imag))

      return

FHKL_BLOCKSIZE=32
mod_fhkl_sorted ="""
__device__ __constant__ double reciprocal_metrical_matrix[6];
__device__ __constant__ double gaussians[%%d];
__device__ __constant__ double symmetry[%%d];

__device__ double get_d_star_sq(double const& H,double const& K,double const& L){
   return
            (H * H) * reciprocal_metrical_matrix[0]
          + (K * K) * reciprocal_metrical_matrix[1]
          + (L * L) * reciprocal_metrical_matrix[2]
          + (2 * H * K) * reciprocal_metrical_matrix[3]
          + (2 * H * L) * reciprocal_metrical_matrix[4]
          + (2 * K * L) * reciprocal_metrical_matrix[5];
}

__device__ double get_H_dot_X(double* const& site, unsigned int const& i,
                              double const& H,double const& K,double const& L,
                              unsigned int const & isym){
  # define SYM_STRIDE %%d
  // case for space group P1
  // return H * site[0] + K * site[1] + L * site[2];

  // but in general,
     return site[3*i+0] * (H * symmetry[SYM_STRIDE*isym]+
                       K * symmetry[SYM_STRIDE*isym+3]+
                       L * symmetry[SYM_STRIDE*isym+6])
          + site[3*i+1] * (H * symmetry[SYM_STRIDE*isym+1]+
                       K * symmetry[SYM_STRIDE*isym+4]+
                       L * symmetry[SYM_STRIDE*isym+7])
          + site[3*i+2] * (H * symmetry[SYM_STRIDE*isym+2]+
                       K * symmetry[SYM_STRIDE*isym+5]+
                       L * symmetry[SYM_STRIDE*isym+8])
          + H * symmetry[SYM_STRIDE*isym+9]
          + K * symmetry[SYM_STRIDE*isym+10]
          + L * symmetry[SYM_STRIDE*isym+11]
    ;
}

__global__ void CUDA_fhkl(double *fhkl_real,double *fhkl_imag,
                        const double *xyzsites, const double *weight,
                        const double *u_iso,
                        const unsigned int scattering_type,
                        const unsigned int index_begin,
                        const unsigned int index_end,
                        const double *hkl)
{
   #define BLOCKSIZE %d
   #define G_STRIDE %%d
   #define USE_DEBYE_WALLER %%d

   const unsigned long ix=threadIdx.x+blockDim.x*blockIdx.x;

   /*later figure out how to make this shared*/
   const double twopi= 8.* atan(1.);
   #if USE_DEBYE_WALLER
     const double eight_pi_sq = 2.* twopi * twopi;
     __shared__ double b_iso[BLOCKSIZE];
   #endif
   /* can the threads coalesce with a stride of 3? Make into __shared__ and See Programming Guide Fig. G-2 */
   const double h = hkl[3*ix];
   const double k = hkl[3*ix+1];
   const double l = hkl[3*ix+2];
   double fr=0,fi=0;

   __shared__ double xyz[3*BLOCKSIZE];
   __shared__ double wght[BLOCKSIZE];

   long at = BLOCKSIZE*index_begin/BLOCKSIZE;
   double x_sq = get_d_star_sq(h,k,l) / 4.;

   /* get unique form factors at d_star_sq.
      all scatterers are the same for any particular kernel invocation; hence only one form factor
    */
   double form_factor = gaussians[G_STRIDE*scattering_type];// constant term
   for (unsigned int term=1; term < G_STRIDE; term+=2){
     form_factor += gaussians[G_STRIDE*scattering_type + term] *
                    exp(-gaussians[G_STRIDE*scattering_type + term + 1] * x_sq);
   }

   for (;at<index_end;at+=BLOCKSIZE) {
      xyz[threadIdx.x]=xyzsites[3*at+threadIdx.x];
      xyz[BLOCKSIZE+threadIdx.x]=xyzsites[3*at+threadIdx.x+BLOCKSIZE];
      xyz[2*BLOCKSIZE+threadIdx.x]=xyzsites[3*at+threadIdx.x+2*BLOCKSIZE];
      wght[threadIdx.x]=weight[at+threadIdx.x];
      #if USE_DEBYE_WALLER
        b_iso[threadIdx.x] = eight_pi_sq * u_iso[at+threadIdx.x];
      #endif
      __syncthreads();

      for(unsigned int i=0;i<BLOCKSIZE;i++) {
        double form_factor_mask = form_factor;
        double debye_waller_factor = 1.;
        #if USE_DEBYE_WALLER
          debye_waller_factor *= exp( -b_iso[i] * x_sq );
        #endif
        if (at+i < index_begin || at+i >= index_end) {
          form_factor_mask = 0.;
          debye_waller_factor = 0.;
        }

        double s,c,s_imag=0.,c_real=0.;
        for(unsigned int isym=0; isym < %%d; isym++){
          /*do not use the faster, but less accurate intrinsic function.*/
          sincos(twopi* get_H_dot_X( xyz, i, h, k, l, isym), &s,&c);
          c_real += c;
          s_imag += s;
        }
        fr += form_factor_mask * wght[i] * c_real * debye_waller_factor;
        fi += form_factor_mask * wght[i] * s_imag * debye_waller_factor;
      }
   }

   fhkl_real[ix] += fr;
   fhkl_imag[ix] += fi;
}
"""%(FHKL_BLOCKSIZE)

class fcalc_container:
  def __init__(self,fcalc):
    self.fcalc = fcalc
  def f_calc(self): return self.fcalc

class direct_summation_simple:
  use_alt_parallel=True
  def __init__(self):
    self.simple=True
  def __eq__(self,other):
    return other=="direct"

class direct_summation_cuda_platform(direct_summation_simple):
  def __init__(self):
    self.pycuda=True
    self.validate_platform_resources()
  def validate_platform_resources(self):
    try:
      import numpy
      self.numpy = numpy
    except:
      raise Sorry("""Module numpy must be installed for parallel processor direct summation.""")

    try:
      import pycuda.autoinit # import dependency
    except:
      raise Sorry("""Module pycuda must be installed for parallel processor direct summation.""")
