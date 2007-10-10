""" Charge flipping algorithm(s) and related data structures

References.

[1] G. Oszl{\'a}nyi and A. S{\"u}t{\H o}.
Ab initio structure solution by charge flipping.
Acta Cryst. A, 60:134--141, 2003.

[2] G. Oszl{\'a}nyi and A. S{\"u}t{\H o}.
Ab initio structure solution by charge flipping.
II. use of weak reflections. Acta Cryst. A, 61:147, 2004.

[3] L. Palatinus and G. Chapuis
SUPERFLIP -- a computer program for the solution of crystal structures
by charge flipping in arbitry dimensions
J. Appl. Cryst., 40:786--790, 2007

[4] M. Shiono and M.M. Woolfson.
Direct-space methods in phase extension and phase determination.
I. low-density elimination. Acta Cryst. A, 48:451-456, 1992.
--> This is a protein paper

[5] H. Takakura, M. Shiono, T.J. Sato, A. Yamamoto, and A.P. Tsai.
Ab initio structure determination of icosahedral zn-mg-ho quasicrystals
by density modification method.
Phys. Rev. Lett., 86:236, 2001
--> This is an elaboration on the method in [4] as well as an application in a
different compartment of crystallography. This is also the method used
in SUPERFLIP circa Sept 2007 to polish the electron density after the charge
flipping method has converged.
"""

from __future__ import division

from libtbx import forward_compatibility
from libtbx import object_oriented_patterns as oop
from libtbx.math_utils import are_equivalent
from libtbx.assert_utils import is_numeric

from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
from cctbx import translation_search

import ab_initio

import math


class _array_extension(oop.injector, miller.array):

  def oszlanyi_suto_phase_transfer(self,
                                   source,
                                   delta_varphi=math.pi/2,
                                   weak_reflection_fraction=0.2,
                                   need_sorting=True):
    """ As per ref. [2] """
    cut = int(weak_reflection_fraction * source.size())
    if need_sorting:
      p = self.sort_permutation(by_value="data", reverse=True)
      target = self.select(p)
      source = source.select(p)
    else:
      target = self
    source_phases = flex.arg(source.data())
    # weak reflections
    phases = source_phases[:cut] + delta_varphi
    moduli = flex.abs(source.data()[:cut])
    # strong ones
    phases.extend(source_phases[cut:])
    moduli.extend(self.data()[cut:])
    return miller.array(self, moduli).phase_transfer(phases)


class _fft_extension(oop.injector, miller.fft_map):
  """ We add those methods to fft_map so that they can be easily reused and
  tested independently of the class charge_flipping_iterator and co. """

  def flipped_fraction_as_delta(self, fraction):
    rho = self.real_map_unpadded().as_1d()
    p = flex.sort_permutation(rho)
    sorted_rho = rho.select(p)
    return sorted_rho[int(fraction * sorted_rho.size())]

  def c_flip(self, delta):
    rho = self.real_map_unpadded().as_1d()
    flipped = rho < delta
    return flex.sum(flex.abs(rho.select(flipped)))

  def c_tot(self):
    return flex.sum(self.real_map())

  def c_tot_over_c_flip_as_delta(self, low, mid, high):
    pass

  def adjusted_delta_based_on_c_tot_over_c_flip(delta):
    r = self.c_tot() / r.c_flip(delta)
    if r < 0.8: return delta * 0.9  # magic numbers from SUPERFLIP
    elif r > 1: return delta * 1.07
    else: return delta



class density_modification_iterator(object):

  def __init__(self, f_obs, f_calc=None, f_000=None):
    assert f_obs.data() is not None
    assert are_equivalent(f_calc is None, f_000 is None)
    assert f_calc is None or (is_numeric(f_calc) and is_numeric(f_000))

    self.original_f_obs = f_obs
    self.f_obs = f_obs.expand_to_p1()\
        .as_non_anomalous_array()\
        .merge_equivalents().array()\
        .discard_sigmas()

    self.symmetry_flags = maptbx.use_space_group_symmetry
    if f_calc is None:
      # Initial f_calc is f_obs with f_000 = 0 and random phases
      f_calc = self.f_obs.randomize_phases()
      f_000 = 0
    self.f_calc = f_calc
    self.f_000 = f_000

  def __iter__(self):
    return self

  def next(self):
    self.compute_electron_density_map()
    self.modify_electron_density()
    self.compute_structure_factors()
    self.transfer_phase_to_f_obs()
    self.f_000 = self._g_000

    return self # iterator-is-its-own-state trick

  def compute_electron_density_map(self):
    self.rho_map = self.f_calc.fft_map(
      f_000=self.f_000,
      symmetry_flags=maptbx.use_space_group_symmetry)

  def compute_structure_factors(self):
    """ This shall compute the structure factors self._g of self.rho_map,
    as well as the 000 component self._g_000, scaling them by the number of
    grid points """
    rho = self.rho_map.real_map()
    scale = 1/self.rho_map.n_grid_points()
    self._g_000 = flex.sum(rho) * scale
    self._g = self.f_obs.structure_factors_from_map(rho, in_place_fft=True)
    self._g *= scale

  def transfer_phase_to_f_obs(self):
    self.f = self.f_obs.phase_transfer(self._g)

  def r1_factor(self):
    return self.f_obs.r1_factor(self._g, assume_index_matching=True)

  def correlation_map_peak_cluster_analysis(self,
                                            return_correlation_map_too=False):
    """ The fast correlation map as per cctbx.translation_search.fast_nv1995
    is computed and its peaks studied """
    f_obs = self.original_f_obs
    f_calc = self.f_calc
    crystal_gridding = f_obs.crystal_gridding(
      symmetry_flags=translation_search.symmetry_flags(
        is_isotropic_search_model=False,
        have_f_part=False))
    correlation_map = translation_search.fast_nv1995(
      gridding=crystal_gridding.n_real(),
      space_group=f_obs.space_group(),
      anomalous_flag=f_obs.anomalous_flag(),
      miller_indices_f_obs=f_obs.indices(),
      f_obs=f_obs.data(),
      f_part=flex.complex_double(), ## no sub-structure is already fixed
      miller_indices_p1_f_calc=f_calc.indices(),
      p1_f_calc=f_calc.data()).target_map()
    ## The correlation map is not a miller.fft_map, just a 3D flex.double
    search_parameters = maptbx.peak_search_parameters(
      peak_search_level=1,
      peak_cutoff=0.1,
      interpolate=True,
      min_distance_sym_equiv=1.,
      max_clusters=24)
    result = crystal_gridding.tags().peak_search(
      map=correlation_map,
      parameters=search_parameters)
    if return_correlation_map_too:
      result = (result, correlation_map)
    return result

  def apply_shift(self, t):
    phase_shifts = -2*math.pi*self.f_calc.indices().as_vec3_double().dot(t)
    phase_shifts *= flex.arg(self.f_calc.data())
    self.f_calc.phase_transfer(phase_shifts)


class basic_iterator(density_modification_iterator):
  """ An iterator over the sequence of electron densities and structure
  factors obtained by repeateadly applying the basic charge flipping
  described in ref. [1].

  Notes.

    self.rho is actually V*rho where V is the unit cell volume
  """

  def __init__(self, f_obs, delta, **kwds):
    super(basic_iterator, self).__init__(f_obs, **kwds)
    self.delta = delta
    if delta is None: self.initialise_delta()

  def initialise_delta(self):
    self.delta = self.rho_map.flipped_fraction_as_delta(0.8)

  def modify_electron_density(self):
    """ This shall modify rho in place """
    ab_initio.ext.flip_charges_in_place(self.rho_map.real_map(), self.delta)


class weak_reflection_improved_iterator(basic_iterator):
  """ The variation described in ref. [2] """

  def __init__(self, f_obs, delta=None,
               delta_varphi=math.pi/2,
               weak_reflection_fraction=0.2,
               **kwds):
    super(weak_reflection_improved_iterator,
          self).__init__(f_obs, delta, **kwds)
    self.delta_varphi = delta_varphi
    self.weak_reflection_fraction = weak_reflection_fraction

    # sort f_obs by increasing amplitudes once and for all
    p = self.f_obs.sort_permutation(by_value="data", reverse=True)
    self.f_obs = self.f_obs.select(p)

  def transfer_phase_from_g_tof_obs(self):
    self.f = self.f_obs.oszlanyi_suto_phase_transfer(
      self._g,
      self.delta_varphi,
      self.weak_reflection_fraction,
      need_sorting=False)


class low_density_elimination_iterator(density_modification_iterator):
  """ A method related to charge flipping, especially useful to sharpen
  the map after the convergence.
  C.f. Ref [4].
  """

  def __init__(self, f_obs,
               rho_c=lambda positives: 0.2*flex.mean(positives),
               **kwds):
    super(low_density_elimination_iterator, self).__init__(f_obs, **kwds)
    self.rho_c = rho_c

  def modify_electron_density(self):
    ab_initio.ext.low_density_elimination_in_place_tanaka_et_al_2001(
      self.rho_map.real_map(), self.rho_c())

  def compute_electron_density_map(self):
    self.rho_map = self.f_calc.fft_map(
      symmetry_flags=maptbx.use_space_group_symmetry)
