""" Charge flipping related tools """

from __future__ import generators, division

from cctbx.array_family import flex
from cctbx import miller
import math

class charge_flipping_iterator(object):
  """
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

  Notes.

    self.rho is actually V*rho where V is the unit cell volume
  """

  def __init__(self, f_obs, delta=None):
    assert f_obs.data() is not None
    self.delta = delta

    self.f_obs = f_obs.expand_to_p1()\
                      .as_non_anomalous_array()\
                      .merge_equivalents().array()\
                      .discard_sigmas()

    self.set_rho_to_randomized_f_obs()

    # delta being None means we need to find a starting value
    if self.delta is None: self.initialise_delta()

  def initialise_delta(self):
    self.set_delta_based_on_flipped_fraction()

  def set_rho_to_randomized_f_obs(self):
    # Initial rho is the alter ego of f_obs with f_000 = 0 and random phases
    random_phases = (2*math.pi)*flex.random_double(self.f_obs.size())
    f = self.f_obs.phase_transfer(random_phases)
    self.rho = f.fft_map().real_map_unpadded()

  def __iter__(self):
    return self

  def next(self):
    # flip
    rho_1d_view = self.rho.as_1d()
    flipped_selection = rho_1d_view < self.delta
    flipped = rho_1d_view.select(flipped_selection)
    flipped *= -1
    rho_1d_view.set_selected(flipped_selection, flipped)

    #FFT (it seems the natural place to scale by the number of grid points)
    scale = 1/self.rho.size()
    self.g = self.f_obs.structure_factors_from_map(self.rho)
    self.g *= scale
    self.g_000 = flex.sum(self.rho) * scale

    f = self.transfer_phases_from_g_to_f_obs()

    # inverse FFT
    self.rho = f.fft_map(f_000=self.g_000).real_map_unpadded()

    # iterator-is-its-own-state trick
    return self

  def transfer_phases_from_g_to_f_obs(self):
    return self.f_obs.phase_transfer(self.g)

  def set_delta_based_on_flipped_fraction(self):
    ## possible optimisation here: partial_sort
    rho = self.rho.as_1d()
    p = flex.sort_permutation(rho)
    sorted_rho = rho.select(p)
    self.delta = sorted_rho[int(sorted_rho.size()*0.8)]

  def set_delta_based_on_c_tot_over_c_flip(self):
    rho = self.rho.as_1d()
    c_tot = flex.sum(rho)
    c_flip = flex.sum(flex.abs(rho).select(rho < self.delta))
    r = c_tot/c_flip
    if r < 0.8: self.delta *= 0.9  # magic numbers from SUPERFLIP
    elif r > 1: self.delta *= 1.07


class charge_flipping_iterator_tweaked_for_weak_reflections(
  charge_flipping_iterator):

  def __init__(self, f_obs,
               delta_varphi=math.pi/2,
               weak_reflection_fraction=0.2
               ):
    super(charge_flipping_iterator_tweaked_for_weak_reflections,
          self).__init__(f_obs)

    assert ( None not in (delta_varphi, weak_reflection_fraction)
             or (delta_varphi, weak_reflection_fraction) == (None, None) )
    self.weak_reflection_fraction = weak_reflection_fraction
    self.delta_varphi = delta_varphi

    self.sorting_permutation = flex.sort_permutation(self.f_obs.data())
    self.weak_cut = int(self.f_obs.size()*self.weak_reflection_fraction)

  def transfer_phases_from_g_to_f_obs(self):
    phases = self.g.phases().data()
    p = self.sorting_permutation
    i = self.weak_cut

    phi_weaks = phases[:i] + self.delta_varphi
    phases.set_selected(p[:i], phi_weaks)

    moduli = flex.double(self.g.size())
    moduli.set_selected(p[:i], self.g.amplitudes().data()[:i])
    moduli.set_selected(p[i:],          self.f_obs.data()[i:])

    return miller.array(self.f_obs, moduli).phase_transfer(phases)
