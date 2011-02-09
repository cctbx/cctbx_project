import copy,math
from scitbx.array_family import flex
from cctbx.xray import ext
from libtbx.utils import Sorry

# alternative implementations of direct summation using parallel architecture
# example usage is given in cctbx/regression/tst_direct_scaling

class pprocess:
  """two things remain for the first version
  2.  limits on # atoms, #hkl etc.
  """

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

  def prepare_data_arrays_for_cuda(self,verbose=False):

    """ The number of miller indices and atoms each must be an exact
        multiple of the BLOCKSIZE (32), so zero-padding is employed"""

    # Transfer data from flex arrays to numpy & zero-pad
    # Marshal the miller indices, fractional coordinates, and occupancies
    # into numpy arrays for use in CUDA

    # Miller indices
    flat_mix = self.miller_indices.as_vec3_double().as_double().as_numpy_array()
    if self.miller_indices.size()%FHKL_BLOCKSIZE > 0:
      newsize = 3*FHKL_BLOCKSIZE* (1+(self.miller_indices.size()/FHKL_BLOCKSIZE))
      if verbose: print "Miller indices: resetting %s array from size %d to %d"%(
                        flat_mix.dtype,flat_mix.shape[0],newsize)
      flat_mix.resize(( newsize,))

    self.n_flat_hkl = flat_mix.shape[0]/3

    # Scattering sites
    sites = self.scatterers.extract_sites()
    n_sites = sites.size()
    flat_sites = sites.as_double().as_numpy_array()
    #print flat_sites.dtype,flat_sites.shape

    # Occupancies
    occupancies = self.scatterers.extract_occupancies().as_numpy_array()
    #print occupancies.dtype,occupancies.shape
    if n_sites%FHKL_BLOCKSIZE > 0:
      newsize = FHKL_BLOCKSIZE* (1+(n_sites/FHKL_BLOCKSIZE))
      if verbose:
        print "Scatterer xyzs: resetting %s array from size %d to %d"%(
              flat_sites.dtype,flat_sites.shape[0],3*newsize)
        print "Occupancies   : resetting %s array from size %d to %d"%(
              occupancies.dtype,occupancies.shape[0],newsize)
        print
      flat_sites.resize(( 3*newsize,))
      # zero-padding for occupancies guarantees no effect from padded sites
      occupancies.resize(( newsize,))

    self.n_sites = n_sites
    self.flat_mix = flat_mix
    self.flat_sites = flat_sites
    self.occupancies = occupancies

  def prepare_registry_symmetry_cell(self,numpy):
    # Marshall the scattering types, form factors, symmetry elements,
    # and metrical matrix into numpy arrays for use in CUDA.

    # Scattering sites registry
    uniqueix = self.registry.unique_indices(self.scatterers
      ).as_numpy_array().astype(numpy.uint32)

    if len(uniqueix)%FHKL_BLOCKSIZE > 0:
      uniqueix.resize(FHKL_BLOCKSIZE* (1+(len(uniqueix)/FHKL_BLOCKSIZE)),)
    assert len(uniqueix) == len(self.occupancies)
    self.uniqueix = uniqueix

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

    # Unit cell dimensions --> metrical matrix --> establishes dstar for given hkl
    self.reciprocal_metrical_matrix = numpy.array(
      self.unit_cell.reciprocal_metrical_matrix() )

  def validate_the_inputs(self,manager,cuda):

    # Assess limits due to the current implementation:
    # no support for space group symmetry, B factors, ADP's, f' or f"

    PyNX_implementation = "Fhkl = Sum(atoms) (occ * exp(2*pi*i(H*X)))"
    # reference http://pynx.sourceforge.net (http://arxiv.org/abs/1010.2641v1)

    Present_implementation = "Fhkl = Sum(atoms) (f0(d*) * occ * exp(2*pi*i(H*X)))"

    Future_implementation = """Fhkl = Sum(atoms(X)) ((f0(d*)+f'+i*f") * occ *
        Sum(symops(SR,ST)) exp(2*pi*i(H*SR*X+H*ST) + DW-factor)
      )"""
    try:
      assert manager._miller_set.space_group_info().type().lookup_symbol()=="P 1"
    except:
      raise Sorry("As presently implemented parallel processor direct summation only supports space group P 1.""")

    try:
      assert self.scatterers.extract_u_iso().count(0.0)==len(self.scatterers)
    except:
      raise Sorry("As presently implemented parallel processor direct summation doesn't support non-zero B factors""")

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

      self.prepare_data_arrays_for_cuda()

      self.prepare_registry_symmetry_cell(algorithm.numpy)

      fhkl_real = algorithm.numpy.zeros((self.n_flat_hkl,),algorithm.numpy.float64)
      fhkl_imag = algorithm.numpy.zeros((self.n_flat_hkl,),algorithm.numpy.float64)

      assert cuda.Device.count() >= 1

      device = cuda.Device(0)
      WARPSIZE=device.get_attribute(cuda.device_attribute.WARP_SIZE) # 32
      MULTIPROCESSOR_COUNT=device.get_attribute(cuda.device_attribute.MULTIPROCESSOR_COUNT)

      mod = SourceModule(mod_fhkl_str%(self.gaussians.shape[0],self.g_stride),
                         options=["-use_fast_math"])

      r_m_m_address = mod.get_global("reciprocal_metrical_matrix")[0]
      cuda.memcpy_htod(r_m_m_address, self.reciprocal_metrical_matrix)

      gaussian_address = mod.get_global("gaussians")[0]
      cuda.memcpy_htod(gaussian_address, self.gaussians)

      CUDA_fhkl = mod.get_function("CUDA_fhkl")

      CUDA_fhkl(cuda.InOut(fhkl_real),
                 cuda.InOut(fhkl_imag),
                 cuda.In(self.flat_sites),
                 cuda.In(self.occupancies),
                 cuda.In(self.uniqueix),
                 algorithm.numpy.int32(self.n_sites),
                 cuda.In(self.flat_mix),
                 block=(FHKL_BLOCKSIZE,1,1),
                 grid=((self.n_flat_hkl/FHKL_BLOCKSIZE,1)))

      flex_fhkl_real = flex.double(fhkl_real[0:len(self.miller_indices)])
      flex_fhkl_imag = flex.double(fhkl_imag[0:len(self.miller_indices)])
      flex_complex = flex.complex_double(flex_fhkl_real,flex_fhkl_imag)

      instance._results = fcalc_container(flex_complex)

      return

FHKL_BLOCKSIZE=32
mod_fhkl_str ="""
__device__ __constant__ double reciprocal_metrical_matrix[6];
__device__ __constant__ double gaussians[%%d];

__device__ double get_d_star_sq(double const& H,double const& K,double const& L){
   return
            (H * H) * reciprocal_metrical_matrix[0]
          + (K * K) * reciprocal_metrical_matrix[1]
          + (L * L) * reciprocal_metrical_matrix[2]
          + (2 * H * K) * reciprocal_metrical_matrix[3]
          + (2 * H * L) * reciprocal_metrical_matrix[4]
          + (2 * K * L) * reciprocal_metrical_matrix[5];
}

__global__ void CUDA_fhkl(double *fhkl_real,double *fhkl_imag,
                        const double *xyzsites, const double *occupancy,
                        const unsigned int *unique_scattering_type,
                        const long natoms,
                        const double *hkl)
{
   #define BLOCKSIZE %d
   #define G_STRIDE %%d

   const unsigned long ix=threadIdx.x+blockDim.x*blockIdx.x;

   /*later figure out how to make this shared*/
   const double twopi= 8.* atan(1.);
   /* can the threads coalesce with a stride of 3? Make into __shared__ and See Programming Guide Fig. G-2 */
   const double h = hkl[3*ix];
   const double k = hkl[3*ix+1];
   const double l = hkl[3*ix+2];
   double fr=0,fi=0;

   __shared__ double xyz[3*BLOCKSIZE];
   __shared__ float occ[BLOCKSIZE];
   __shared__ unsigned int unique_type[BLOCKSIZE];
   long at=0;
   double x_sq = get_d_star_sq(h,k,l) / 4.;

   for (;at<natoms;at+=BLOCKSIZE) {
      xyz[threadIdx.x]=xyzsites[3*at+threadIdx.x];
      xyz[BLOCKSIZE+threadIdx.x]=xyzsites[3*at+threadIdx.x+BLOCKSIZE];
      xyz[2*BLOCKSIZE+threadIdx.x]=xyzsites[3*at+threadIdx.x+2*BLOCKSIZE];
      occ[threadIdx.x]=occupancy[at+threadIdx.x];
      unique_type[threadIdx.x]=unique_scattering_type[at+threadIdx.x];
      __syncthreads();

      for(unsigned int i=0;i<BLOCKSIZE;i++) {

        /* take a huge hit here; until such time that we pre-sort the scattering
           types, we'll have to recalculate the form_factor(x_sq) for each
           scatterer.
           Future version: presort the scatterers by atom type */

        /* get unique form factors at d_star_sq */
        double form_factor = gaussians[G_STRIDE*unique_type[i]];// constant term
        for (unsigned int term=1; term < G_STRIDE; term+=2){
          form_factor += gaussians[G_STRIDE*unique_type[i] + term] *
                     exp(-gaussians[G_STRIDE*unique_type[i] + term + 1] * x_sq);
        }
        if (at+i>=natoms) form_factor = 0.;

        double s,c;
        /*do not use the faster, but less accurate intrinsic function.*/
        sincos(twopi* (h*xyz[3*i] + k*xyz[3*i+1] + l*xyz[3*i+2]) , &s,&c);
        fr += form_factor * occ[i]*c;
        fi += form_factor * occ[i]*s;
      }
   }

   fhkl_real[ix]+=fr;
   fhkl_imag[ix]+=fi;
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
