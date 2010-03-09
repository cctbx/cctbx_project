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

from __future__ import division, generators

from libtbx import forward_compatibility
from libtbx import object_oriented_patterns as oop
from libtbx.math_utils import are_equivalent
from libtbx.assert_utils import is_numeric
from libtbx import adopt_init_args, adopt_optional_init_args

from scitbx import matrix as mat

from cctbx.array_family import flex
from cctbx import crystal
from cctbx import sgtbx
from cctbx import miller
from cctbx import maptbx
from cctbx import translation_search

from smtbx import ab_initio

import scitbx.math

import itertools
import sys
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
  tested independently of the charge flipping iterators. """

  def flipped_fraction_as_delta(self, fraction):
    rho = self.real_map_unpadded().as_1d()
    p = flex.sort_permutation(rho)
    sorted_rho = rho.select(p)
    return sorted_rho[int(fraction * sorted_rho.size())]
  flipped_fraction_as_delta = oop.memoize_method(flipped_fraction_as_delta)

  def c_flip(self, delta):
    rho = self.real_map_unpadded().as_1d()
    return flex.sum(flex.abs(rho.select(rho < delta)))
  c_flip = oop.memoize_method(c_flip)

  def c_tot(self):
    return flex.sum(self.real_map())
  c_tot = oop.memoize_method(c_tot)

  def skewness(self):
    return maptbx.more_statistics(self.real_map()).skewness()
  skewness = oop.memoize_method(skewness)

  def sigma(self):
    return maptbx.statistics(self.real_map()).sigma()
  sigma = oop.memoize_method(sigma)


class density_modification_iterator(object):
  """ Skeleton for any method which, like charge flipping, does cycles like

      rho --|1|--> rho' --|Fourier analysis|--> g --|2|--> f
       ^                                                   |
       |----------------|Fourier synthesis|----------------|

    where the transformation (1) and (2) are specific to each method.

    Synopsis:
      flipping = heir_of_density_modification_iterator(...)
      flipping.start(f_obs, initial_phases)
      flipping.next() # 1st cycle
      flipping.next() # 2nd cycle
      ....
  """

  min_cc_peak_height = 0.8

  def __init__(self, **kwds):
    adopt_optional_init_args(self, kwds)

  def start(self, f_obs, phases, f_000=0):
    self.f_obs = f_obs
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell=self.f_obs.unit_cell(),
      space_group_info=sgtbx.space_group_info('P1'),
      d_min=self.f_obs.d_min(),
      resolution_factor=1/2,
      symmetry_flags=maptbx.use_space_group_symmetry)

    self.fft_scale = (self.f_obs.crystal_symmetry().unit_cell().volume()
                      / self.crystal_gridding.n_grid_points())
    self.f_calc = self.f_obs.phase_transfer(phases)
    self.f_000 = f_000
    self.compute_electron_density_map()

  def __iter__(self):
    return self

  def next(self):
    """ perform one cycle and return itself """
    self.modify_electron_density()
    self.compute_structure_factors()
    self.transfer_phase_to_f_obs()
    self.f_000 = self._g_000
    self.compute_electron_density_map()
    return self # iterator-is-its-own-state trick

  def compute_electron_density_map(self):
    """ Compute the electron density from the structure factors self.f_calc
    and the 000 component self.f_000, scaling by the unit cell volume """
    self.rho_map = miller.fft_map(self.crystal_gridding,
                                  self.f_calc,
                                  self.f_000)
    self.rho_map.apply_volume_scaling()

  def compute_structure_factors(self):
    """ Compute the structure factors self._g of self.rho_map,
    as well as the 000 component self._g_000, scaling them by the number of
    grid points """
    rho = self.rho_map.real_map()
    self._g_000 = flex.sum(rho) * self.fft_scale
    self._g = self.f_obs.structure_factors_from_map(rho, in_place_fft=True)
    self._g *= self.fft_scale

  def transfer_phase_to_f_obs(self):
    self.f_calc = self.f_obs.phase_transfer(self._g)

  def r1_factor(self):
    return self.f_obs.r1_factor(self._g, assume_index_matching=True)


class basic_iterator(density_modification_iterator):
  """ An iterator over the sequence of electron densities and structure
  factors obtained by repeateadly applying the basic charge flipping
  described in ref. [1].
  """

  def __init__(self, delta=None, **kwds):
    super(basic_iterator, self).__init__(**kwds)
    self.delta = delta

  def c_tot_over_c_flip(self):
    return self.rho_map.c_tot()/self.rho_map.c_flip(self.delta)

  def modify_electron_density(self):
    """ This shall modify rho in place """
    ab_initio.ext.flip_charges_in_place(self.rho_map.real_map(), self.delta)


class weak_reflection_improved_iterator(basic_iterator):
  """ The variation described in ref. [2] """

  def __init__(self, delta=None,
               delta_varphi=math.pi/2,
               weak_reflection_fraction=0.2,
               **kwds):
    super(weak_reflection_improved_iterator,
          self).__init__(delta, **kwds)
    self.delta_varphi = delta_varphi
    self.weak_reflection_fraction = weak_reflection_fraction

  def start(self, f_obs, phases, f_000=0):
    """ sort f_obs by increasing amplitudes once and for all """
    super(weak_reflection_improved_iterator, self).start(f_obs, phases, f_000)
    p = self.f_obs.sort_permutation(by_value="data", reverse=True)
    self.f_obs = self.f_obs.select(p)

  def transfer_phase_to_f_obs(self):
    self.f_calc = self.f_obs.oszlanyi_suto_phase_transfer(
      self._g,
      self.delta_varphi,
      self.weak_reflection_fraction,
      need_sorting=False)


class low_density_elimination_iterator(density_modification_iterator):
  """ A method related to charge flipping.
  C.f. Ref [4].
  """

  def __init__(self, constant_rho_c=None, **kwds):
    super(low_density_elimination_iterator, self).__init__(**kwds)
    self.constant_rho_c = constant_rho_c

  def modify_electron_density(self):
    ab_initio.ext.low_density_elimination_in_place_tanaka_et_al_2001(
      self.rho_map.real_map(), self.rho_c())

  def rho_c(self):
    if self.constant_rho_c is not None:
      return self.constant_rho_c
    else:
      return self.shiono_woolfson_rho_c()

  def shiono_woolfson_rho_c(self):
    """ The rho_c suggested in Ref [4] """
    rho = self.rho_map.real_map_unpadded().as_1d()
    return 0.2*flex.mean(rho.select(rho >0))


def f_calc_symmetrisations(f_obs, f_calc_in_p1, min_cc_peak_height):
  # The fast correlation map as per cctbx.translation_search.fast_nv1995
  # is computed and its peaks studied.
  # Inspiration from phenix.substructure.hyss for the parameters tuning.
  crystal_gridding = f_obs.crystal_gridding(
    symmetry_flags=translation_search.symmetry_flags(
      is_isotropic_search_model=False,
      have_f_part=False),
    resolution_factor=1/3
  )
  correlation_map = translation_search.fast_nv1995(
    gridding=crystal_gridding.n_real(),
    space_group=f_obs.space_group(),
    anomalous_flag=f_obs.anomalous_flag(),
    miller_indices_f_obs=f_obs.indices(),
    f_obs=f_obs.data(),
    f_part=flex.complex_double(), ## no sub-structure is already fixed
    miller_indices_p1_f_calc=f_calc_in_p1.indices(),
    p1_f_calc=f_calc_in_p1.data()).target_map()

  if 0:
    from crys3d import wx_map_viewer
    wx_map_viewer.display(title="CC map",
                          raw_map=correlation_map,
                          unit_cell=f_calc_in_p1.unit_cell())

  search_parameters = maptbx.peak_search_parameters(
    peak_search_level=1,
    peak_cutoff=0.5,
    interpolate=True,
    min_distance_sym_equiv=1e-6,
    general_positions_only=False,
    min_cross_distance=f_obs.d_min()/2)
  ## The correlation map is not a miller.fft_map, just a 3D flex.double
  correlation_map_peaks = crystal_gridding.tags().peak_search(
    map=correlation_map,
    parameters=search_parameters)
  # iterate over the strong peak; for each, shift and symmetrised f_calc
  for peak in correlation_map_peaks:
    if peak.height < min_cc_peak_height: break
    shift = mat.col(peak.site)
    multiplier = int(1e4)
    shift_1 = shift*multiplier
    cb_op = sgtbx.change_of_basis_op(
      sgtbx.rt_mx(sgtbx.tr_vec(shift_1.as_int(), multiplier)))
    symm_f_calc = (f_calc_in_p1
      .change_basis(cb_op)
      .customized_copy(
        space_group_info=f_obs.space_group_info()))
    merging = symm_f_calc.merge_equivalents()
    symm_f_calc = merging.array()
    yield symm_f_calc, shift, peak.height


class solving_iterator(object):

  initial_phases_type = "random"
  delta_guessing_method = "sigma"
  delta_over_sigma = 1.1
  max_delta_guessing_iterations = 10
  map_sigma_stability_threshold = 0.01
  initial_flipped_fraction=0.8
  yield_during_delta_guessing = False
  max_solving_iterations = 500
  max_attempts_to_get_phase_transition = 5
  max_attempts_to_get_sharp_correlation_map = 5
  yield_solving_interval = 10
  polishing_iterations = 5
  min_cc_peak_height = 0.9

  def __init__(self, flipping_iterator, f_obs, **kwds):
    self.flipping_iterator = flipping_iterator
    adopt_optional_init_args(self, kwds)
    self.attempts = []
    self.f_calc_solutions = []
    self.had_phase_transition = False
    self.max_attempts_exceeded = False
    self.state = self.starting = {
      "random": self._starting_with_random_phases,
      }[self.initial_phases_type](f_obs)
    self.guessing_delta = {
      "sigma": self._guessing_delta_with_map_sigma,
      "c_tot_over_c_flip": self._guessing_delta_with_c_tot_over_c_flip,
      }[self.delta_guessing_method]()
    self.solving = self._solving()
    self.polishing = self._polishing()
    self.evaluating = self._evaluating(f_obs)
    self.finished = self._finished()

  def __iter__(self):
    """ Note: a loop for flipping in solving_iterator_obj: that is
    interrupted by break will reliably result in a call
    solving_iterator_obj.clean_up() in Python 2.5+ while the code should
    still run on earlier versions of Python but without the clean-up. """
    while 1:
      try: state = self.state.next()
      except StopIteration: break
      try: yield self.flipping_iterator
      except GeneratorExit: break
      self.state = state
    self.clean_up()

  def clean_up(self):
    """ The generator-based state machine pattern used to implement this
    class creates cycles for each generator:
       self.polishing.gi_frame.f_locals['self'] is self == True
    for example.
    Thus reference counting does not have it collected,
    and self.flipping_iterator is not collected either.
    The latter holds large objects (a fft_map and a miller.array), which
    results in the memory being used to creep up each time a charge
    flipping run is done.
    Note: using a weak reference for solving.flipping_iterator would not work
    because that object is also owned by several of the generators' frame
    mentionned above.

    Thus we delete the generators after the run has finished, therefore
    breaking the cycle.
    """
    del self.state
    del self.starting
    del self.guessing_delta
    del self.solving
    del self.polishing
    del self.evaluating
    del self.finished

  def _starting_with_random_phases(self, f_obs):
    f_obs = f_obs.eliminate_sys_absent() \
                 .expand_to_p1() \
                 .as_non_anomalous_array() \
                 .merge_equivalents().array() \
                 .discard_sigmas()
    while 1:
      random_phases = (2*math.pi)*flex.random_double(f_obs.size())
      self.flipping_iterator.start(f_obs, random_phases)
      yield self.guessing_delta

  def _finished(self):
    if not self.max_attempts_exceeded:
      self.had_phase_transition = True
    yield self.finished

  def _guessing_delta_with_c_tot_over_c_flip(self):
    flipping = self.flipping_iterator
    delta_needs_initialisation = True
    while 1:
      self.f_calc_solutions = []
      if delta_needs_initialisation:
        flipping.delta = flipping.rho_map.flipped_fraction_as_delta(
                                                self.initial_flipped_fraction)
        delta_needs_initialisation = False
      for foo in itertools.islice(flipping,
                                  self.max_delta_guessing_iterations):
        pass
      r = flipping.c_tot_over_c_flip()
      # magic numbers from SUPERFLIP
      low, high = 0.8, 1.
      if low <= r <= high:
        yield self.solving
        flipping.restart()
        delta_needs_initialisation = True
      else:
        if self.yield_during_delta_guessing:
          yield self.guessing_delta
        if r < low:
          flipping.delta *= 0.9
        elif r > high:
          flipping.delta *= 1.07

  def _guessing_delta_with_map_sigma(self):
    while 1:
      self.f_calc_solutions = []
      sigmas = flex.double()
      for i in xrange(self.max_delta_guessing_iterations):
        sigma = self.flipping_iterator.rho_map.sigma()
        sigmas.append(sigma)
        self.flipping_iterator.delta = self.delta_over_sigma * sigma
        if len(sigmas) == 1:
          self.flipping_iterator.next()
          continue
        sigma_tail_stats = scitbx.math.basic_statistics(sigmas[-5:])
        if (abs(sigma_tail_stats.bias_corrected_standard_deviation
                /sigma_tail_stats.mean) < self.map_sigma_stability_threshold):
          break
        if self.yield_during_delta_guessing: yield self.guessing_delta
        self.flipping_iterator.next()
      yield self.solving

  def _solving(self):
    while 1:
      i_attempt = 0
      while i_attempt < self.max_attempts_to_get_phase_transition:
        i_attempt += 1
        if i_attempt > 2:
          self.max_solving_iterations *= 1.5
        self.skewness_evolution = observable_evolution()
        for n, flipping in enumerate(
          itertools.islice(self.flipping_iterator,
                           0, self.max_solving_iterations)):
          self.iteration_index = n
          if n % self.yield_solving_interval == 0:
            yield self.solving
          self.skewness_evolution.append(flipping.rho_map.skewness())
          #if flipping.rho_map.skewness() < 3: continue
          if self.skewness_evolution.had_phase_transition():
            self.attempts.append(n)
            yield self.polishing
            break
        else:
          if i_attempt != self.max_attempts_to_get_phase_transition:
            yield self.starting
      self.max_attempts_exceeded = True
      yield self.finished

  def _polishing(self):
    while 1:
      low_density_elimination = low_density_elimination_iterator(
        constant_rho_c=self.flipping_iterator.delta)
      low_density_elimination.start(f_obs=self.flipping_iterator.f_obs,
                                    phases=self.flipping_iterator.f_calc,
                                    f_000=0)
      for i in xrange(self.polishing_iterations):
        low_density_elimination.next()
      self.flipping_iterator.f_calc = low_density_elimination.f_calc
      self.flipping_iterator.f_000 = low_density_elimination.f_000
      yield self.evaluating

  def _evaluating(self, f_obs):
    while 1:
      attempts = 0
      while attempts < self.max_attempts_to_get_sharp_correlation_map:
        attempts += 1
        self.f_calc_solutions = []
        for f_calc, shift, cc_peak_height\
            in f_calc_symmetrisations(f_obs,
                                      self.flipping_iterator.f_calc,
                                      self.min_cc_peak_height):
          if (not self.f_calc_solutions
              and cc_peak_height < self.min_cc_peak_height):
            yield self.starting
            break
          self.f_calc_solutions.append((f_calc, shift, cc_peak_height))
        else:
          yield self.finished
      self.max_attempts_exceeded = True

def loop(solving, verbose=True, stdout=sys.stdout):
  previous_state = None
  for flipping in solving:
    if solving.state is solving.guessing_delta:
      # Guessing a value of delta leading to subsequent good convergence
      if verbose:
        if previous_state is solving.solving:
          print "** Restarting (no phase transition) **"
        elif previous_state is solving.evaluating:
          print "** Restarting (no sharp correlation map) **"
      if verbose == "highly":
        if previous_state is not solving.guessing_delta:
          print "Guessing delta..."
          print ("%10s | %10s | %10s | %10s | %10s | %10s | %10s"
                 % ('delta', 'delta/sig', 'R', 'F000',
                    'c_tot', 'c_flip', 'c_tot/c_flip'))
          print "-"*90
        rho = flipping.rho_map
        c_tot = rho.c_tot()
        c_flip = rho.c_flip(flipping.delta)
        # to compare with superflip output
        c_tot *= flipping.fft_scale; c_flip *= flipping.fft_scale
        print "%10.4f | %10.4f | %10.3f | %10.3f | %10.1f | %10.1f | %10.2f"\
              % (flipping.delta, flipping.delta/rho.sigma(),
                 flipping.r1_factor(), flipping.f_000,
                 c_tot, c_flip, c_tot/c_flip)

    elif solving.state is solving.solving:
      # main charge flipping loop to solve the structure
      if verbose=="highly":
        if previous_state is not solving.solving:
          print
          print "Solving..."
          print "with delta=%.4f" % flipping.delta
          print
          print "%5s | %10s | %10s" % ('#', 'F000', 'skewness')
          print '-'*33
        print "%5i | %10.1f | %10.3f" % (solving.iteration_index,
                                         flipping.f_000,
                                         flipping.rho_map.skewness())

    elif solving.state is solving.polishing:
      if verbose == 'highly':
        print
        print "Polishing"
    elif solving.state is solving.finished:
      if solving.max_attempts_exceeded:
        print
        print "** Maximum number of attempts exceeded: it won't solve!"
      break
    previous_state = solving.state


class observable_evolution(object):

  smoothing_coefficient = 0.25
  increasing = True
  noise_level_before = 0.3
  noise_level_after = 0.2

  def __init__(self, **kwds):
    adopt_optional_init_args(self, kwds)
    self.values = flex.double()
    self.raw_values = flex.double()
    self.differences = flex.double()

  def append(self, x):
    self.raw_values.append(x)
    if len(self.values) > 1:
      a = self.smoothing_coefficient
      y0 = self.values[-1]
      y1 = y0 + a*(x - y0)
      self.values.append(y1)
      delta = y1 - y0
      if not self.increasing: delta = -delta
      if len(self.differences) > 1:
        eta0 = self.differences[-1]
        eta1 = eta0 + a*(delta - eta0)
        self.differences.append(eta1)
      else:
        self.differences.append(delta)
    else:
      self.values.append(x)

  def had_phase_transition(self):
    try:
      if len(self.differences) < 5: return False
      i_max = flex.max_index(self.differences)
      max_difference = self.differences[i_max]
      if len(self.differences) - i_max < 10: return False
      i_after  = flex.first_index(self.differences[i_max:] > 0, False)
      if i_after is None: return False
      i_after += i_max
      if len(self.differences) - i_after < 20:
        return False
      i_before = flex.last_index(self.differences[:i_max] > 0, False)
      if i_before is None:
        # observable kept increasing from start to candidate transition
        # find inflexion point
        i_before = flex.last_index(
          self.differences[1:i_max] > self.differences[:i_max-1],
          False)
        if i_before is None: return False
      before = self.differences[2:i_before] # skip first iterations
      min_before = max(flex.min(before), 0)
      after  = self.differences[(i_after - len(self.differences))//2:]
      noisy_before = (
        flex.abs(before - min_before)
        > self.noise_level_before * (max_difference - min_before))
      noisy_after  = (
        flex.abs(after) > self.noise_level_after  * max_difference)
      if noisy_before.count(True)/len(before) > 0.1:
        return False
      if noisy_after.count(True)/len(after) > 0.1:
        return False
      return True
    except ZeroDivisionError:
      return False
