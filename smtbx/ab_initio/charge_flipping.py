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

"""

from __future__ import division

from libtbx import object_oriented_patterns as oop
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
import math

class _array_extension(oop.extends, miller.array):

  def oszlanyi_suto_phase_transfer(self,
                                   source,
                                   delta_varphi=math.pi/2,
                                   weak_reflection_fraction=0.2,
                                   need_sorting=True):
    """ As per ref. [2] """
    cut = int(weak_reflection_fraction * source.size())
    if need_sorting:
      p = self.sort_permutation(by_value="data")
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


class _fft_extension(oop.extends, miller.fft_map):
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


class basic_charge_flipping_iterator(object):
  """ An iterator over the sequence of electron densities and structure
  factors obtained by repeateadly applying the basic charge flipping
  described in ref. [1].

  Notes.

    self.rho is actually V*rho where V is the unit cell volume
  """

  def __init__(self, f_obs, delta=None, relative_delta=None):
    assert f_obs.data() is not None
    assert delta is None or relative_delta is None
    if relative_delta is not None:
      delta = relative_delta * flex.max(f_obs.data())
    self.delta = delta

    self.f_obs = f_obs.expand_to_p1()\
                      .as_non_anomalous_array()\
                      .merge_equivalents().array()\
                      .discard_sigmas()

    # Initial rho is the alter ego of f_obs with f_000 = 0 and random phases
    f = self.f_obs.randomize_phases()
    self.symmetry_flags = maptbx.use_space_group_symmetry
    self.rho_map = f.fft_map(symmetry_flags=self.symmetry_flags)

    # delta being None means we need to find a starting value
    if self.delta is None: self.initialise_delta()

  def initialise_delta(self):
    self.delta = self.rho_map.real_map_unpadded().flipped_fraction_as_delta()

  def __iter__(self):
    return self

  def next(self):
    rho = self.rho_map.real_map_unpadded()

   # flip
    rho_1d_view = rho.as_1d()
    flipped_selection = rho_1d_view < self.delta
    flipped = rho_1d_view.select(flipped_selection)
    flipped *= -1
    rho_1d_view.set_selected(flipped_selection, flipped)

    #FFT (it seems the natural place to scale by the number of grid points)
    scale = 1/rho.size()
    self.g = self.f_obs.structure_factors_from_map(rho)
    self.g *= scale
    self.g_000 = flex.sum(rho) * scale

    # Transfer the phases of g onto f_obs
    f = self.transfer_phase_from_g_to_f_obs()

    # inverse FFT
    self.rho_map = f.fft_map(f_000=self.g_000,
                             symmetry_flags=maptbx.use_space_group_symmetry)

    # iterator-is-its-own-state trick
    return self

  def transfer_phase_from_g_to_f_obs(self):
    return self.f_obs.phase_transfer(self.g)

  def r1_factor(self):
    return self.f_obs.r1_factor(self.g, assume_index_matching=True)


class weak_reflection_improved_charge_flipping_iterator(
  basic_charge_flipping_iterator):
  """ The variation described in ref. [2] """

  def __init__(self, f_obs, delta=None,
                     delta_varphi=math.pi/2,
                     weak_reflection_fraction=0.2):
    super(weak_reflection_improved_charge_flipping_iterator,
          self).__init__(f_obs, delta)
    self.delta_varphi = delta_varphi
    self.weak_reflection_fraction = weak_reflection_fraction

    # sort f_obs by increasing amplitudes once and for all
    p = self.f_obs.sort_permutation(by_value="data", reverse=True)
    self.f_obs = self.f_obs.select(p)

  def transfer_phase_from_g_to_f_obs(self):
    return self.f_obs.oszlanyi_suto_phase_transfer(
      self.g,
      self.delta_varphi,
      self.weak_reflection_fraction,
      need_sorting=False)
