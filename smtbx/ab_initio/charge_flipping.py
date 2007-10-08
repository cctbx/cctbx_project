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

[4] M.ÊShiono and M.M. Woolfson.
Direct-space methods in phase extension and phase determination.
I. low-density elimination. Acta Cryst. A, 48:451Ð456, 1992.
--> This is a protein paper

[5] H.ÊTakakura, M.ÊShiono, T.J. Sato, A.ÊYamamoto, and A.P. Tsai.
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
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
from cctbx import translation_search
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

  def __init__(self, f_obs):
    assert f_obs.data() is not None
    self.original_f_obs = f_obs
    self.f_obs = f_obs.expand_to_p1()\
        .as_non_anomalous_array()\
        .merge_equivalents().array()\
        .discard_sigmas()

    # Initial rho is the alter ego of f_obs with f_000 = 0 and random phases
    f = self.f_obs.randomize_phases()
    self.symmetry_flags = maptbx.use_space_group_symmetry
    self.rho_map = f.fft_map(symmetry_flags=self.symmetry_flags)

  def __iter__(self):
    return self

  def next(self):
    rho = self.rho_map.real_map_unpadded()
    self.modify_electron_density(rho)
    self.compute_structure_factors(rho)
    f = self.transfer_phase_from_g_to_f_obs()
    self.rho_map = f.fft_map(f_000=self.g_000,
                             symmetry_flags=maptbx.use_space_group_symmetry)
    return self # iterator-is-its-own-state trick

  def compute_structure_factors(self, rho):
    """ This shall compute the structure factors self.g of rho,
    as well as the 000 component self.g_000, scaling them by the number of
    grid points """
    scale = 1/rho.size()
    self.g = self.f_obs.structure_factors_from_map(rho)
    self.g *= scale
    self.g_000 = flex.sum(rho) * scale

  def transfer_phase_from_g_to_f_obs(self):
    return self.f_obs.phase_transfer(self.g)

  def r1_factor(self):
    return self.f_obs.r1_factor(self.g, assume_index_matching=True)

  def correlation_map_peak_cluster_analysis(self):
    """ The fast correlation map as per cctbx.translation_search.fast_nv1995
    is computed and its peaks studied """
    f_obs = self.original_f_obs
    f_calc = self.g
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
      min_distance_sym_equiv=1.5,
      max_clusters=24)
    result = crystal_gridding.tags().peak_search(
      map=correlation_map,
      parameters=search_parameters)
    return result


class basic_iterator(density_modification_iterator):
  """ An iterator over the sequence of electron densities and structure
  factors obtained by repeateadly applying the basic charge flipping
  described in ref. [1].

  Notes.

    self.rho is actually V*rho where V is the unit cell volume
  """

  def __init__(self, f_obs, delta):
    super(basic_iterator, self).__init__(f_obs)
    self.delta = delta
    if delta is None: self.initialise_delta()

  def initialise_delta(self):
    self.delta = self.rho_map.flipped_fraction_as_delta(0.8)

  def modify_electron_density(self, rho):
    """ This shall modify rho in place """
    rho_1d_view = rho.as_1d()
    flipped_selection = rho_1d_view < self.delta
    flipped = rho_1d_view.select(flipped_selection)
    flipped *= -1
    rho_1d_view.set_selected(flipped_selection, flipped)


class weak_reflection_improved_iterator(basic_iterator):
  """ The variation described in ref. [2] """

  def __init__(self, f_obs, delta=None,
               delta_varphi=math.pi/2,
               weak_reflection_fraction=0.2):
    super(weak_reflection_improved_iterator, self).__init__(f_obs, delta)
    self.delta_varphi = delta_varphi
    self.weak_reflection_fraction = weak_reflection_fraction

    # sort f_obs by increasing amplitudes once and for all
    p = self.f_obs.sort_permutation(by_value="data", reverse=True)
    self.f_obs = self.f_obs.select(p)

  def transfer_phase_from_g_tof_obs(self):
    return self.f_obs.oszlanyi_suto_phase_transfer(
      self.g,
      self.delta_varphi,
      self.weak_reflection_fraction,
      need_sorting=False)


class low_density_elimination_iterator(density_modification_iterator):
  """ A method related to charge flipping, especially useful to sharpen
  the map after the convergence.
  C.f. Ref [4].
  """

  def __init__(self, f_obs, average_fraction=0.2):
    super(low_density_elimination_iterator, self).__init__(f_obs)
    self.average_fraction = average_fraction

  def modify_electron_density(self, rho):
    rho_1d_view = rho.as_1d()
    negative_selection = rho_1d_view < 0
    rho_1d_view.set_selected(negative_selection, 0)
    positive_selection = ~negative_selection
    positives = rho_1d_view.select(positive_selection)
    rho_c = flex.mean(positives)
    a = -1/(2*(self.average_fraction*rho_c)**2)
    positives *= 1 - flex.exp(a*positives*positives)
    rho_1d_view.set_selected(positive_selection, positives)
