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
from libtbx import itertbx
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

  def __init__(self, f_obs, f_calc=None, f_000=None,
               starter=miller.array.randomize_phases):
    assert f_obs.data() is not None
    assert are_equivalent(f_calc is None, f_000 is None)
    assert f_calc is None or (is_numeric(f_calc) and is_numeric(f_000))

    self.original_f_obs = f_obs
    self.f_obs = f_obs.eliminate_sys_absent() \
                      .expand_to_p1() \
                      .as_non_anomalous_array() \
                      .merge_equivalents().array() \
                      .discard_sigmas()

    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell=self.f_obs.unit_cell(),
      space_group_info=sgtbx.space_group_info('P1'),
      d_min=self.f_obs.d_min(),
      resolution_factor=1/2,
      symmetry_flags=maptbx.use_space_group_symmetry)

    self.fft_scale = (self.f_obs.crystal_symmetry().unit_cell().volume()
                      / self.crystal_gridding.n_grid_points())
    self.starter = starter
    self.restart(f_calc, f_000)

  def restart(self, f_calc=None, f_000=None):
    assert are_equivalent(f_calc is None, f_000 is None)
    assert f_calc is None or (is_numeric(f_calc) and is_numeric(f_000))
    if f_calc is None:
      f_calc = self.starter(self.f_obs)
      f_000 = 0
    self.f_calc = f_calc
    self.f_000 = f_000
    self.compute_electron_density_map()

  def __iter__(self):
    return self

  def next(self):
    self.modify_electron_density()
    self.compute_structure_factors()
    self.transfer_phase_to_f_obs()
    self.f_000 = self._g_000
    self.compute_electron_density_map()
    return self # iterator-is-its-own-state trick

  def compute_electron_density_map(self):
    self.rho_map = miller.fft_map(self.crystal_gridding,
                                  self.f_calc,
                                  self.f_000)
    self.rho_map.apply_volume_scaling()

  def compute_structure_factors(self):
    """ This shall compute the structure factors self._g of self.rho_map,
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

  def f_calc_symmetrisations(self):
    # The fast correlation map as per cctbx.translation_search.fast_nv1995
    # is computed and its peaks studied.
    # Inspiration from phenix.substructure.hyss for the parameters tuning.
    crystal_gridding = self.original_f_obs.crystal_gridding(
      symmetry_flags=translation_search.symmetry_flags(
        is_isotropic_search_model=False,
        have_f_part=False),
      resolution_factor=1/3
    )
    correlation_map = translation_search.fast_nv1995(
      gridding=crystal_gridding.n_real(),
      space_group=self.original_f_obs.space_group(),
      anomalous_flag=self.original_f_obs.anomalous_flag(),
      miller_indices_f_obs=self.original_f_obs.indices(),
      f_obs=self.original_f_obs.data(),
      f_part=flex.complex_double(), ## no sub-structure is already fixed
      miller_indices_p1_f_calc=self.f_calc.indices(),
      p1_f_calc=self.f_calc.data()).target_map()

    if 0:
      from crys3d import wx_map_viewer
      wx_map_viewer.display(title="CC map",
                            raw_map=correlation_map,
                            unit_cell=self.f_calc.unit_cell())

    search_parameters = maptbx.peak_search_parameters(
      peak_search_level=1,
      peak_cutoff=0.5,
      interpolate=True,
      min_distance_sym_equiv=1e-6,
      general_positions_only=False,
      min_cross_distance=self.original_f_obs.d_min()/2)
    ## The correlation map is not a miller.fft_map, just a 3D flex.double
    correlation_map_peaks = crystal_gridding.tags().peak_search(
      map=correlation_map,
      parameters=search_parameters)
    # iterate over the strong peak; for each, shift and symmetrised f_calc
    for peak in correlation_map_peaks:
      if peak.height < 0.8: break
      shift = mat.col(peak.site)
      multiplier = int(1e4)
      shift_1 = shift*multiplier
      cb_op = sgtbx.change_of_basis_op(
        sgtbx.rt_mx(sgtbx.tr_vec(shift_1.as_int(), multiplier)))
      symm_f_calc = (self.f_calc
        .change_basis(cb_op)
        .customized_copy(
          space_group_info=self.original_f_obs.space_group_info()))
      merging = symm_f_calc.merge_equivalents()
      symm_f_calc = merging.array()
      yield symm_f_calc, shift, peak.height


class basic_iterator(density_modification_iterator):
  """ An iterator over the sequence of electron densities and structure
  factors obtained by repeateadly applying the basic charge flipping
  described in ref. [1].
  """

  def __init__(self, f_obs, delta=None, **kwds):
    super(basic_iterator, self).__init__(f_obs, **kwds)
    self.delta = delta

  def c_tot_over_c_flip(self):
    return self.rho_map.c_tot()/self.rho_map.c_flip(self.delta)

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

  def __init__(self, f_obs, constant_rho_c=None, **kwds):
    super(low_density_elimination_iterator, self).__init__(f_obs, **kwds)
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


class solving_iterator(object):

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

  def __init__(self, flipping_iterator, **kwds):
    self.flipping_iterator = flipping_iterator
    adopt_optional_init_args(self, kwds)
    self.attempts = []
    self.f_calc_solutions = []
    self.max_attempts_exceeded = False
    self.state = self.guessing_delta = {
      "sigma": self._guessing_delta_with_map_sigma,
      "c_tot_over_c_flip": self._guessing_delta_with_c_tot_over_c_flip,
      }[self.delta_guessing_method]()
    self.solving = self._solving()
    self.polishing = self._polishing()
    self.evaluating = self._evaluating()
    self.finished = self._finished()

  def __iter__(self):
    while 1:
      try: state = self.state.next()
      except StopIteration: break
      yield self.flipping_iterator
      self.state = state

  def had_phase_transition(self):
    return self.state == self.finished and not self.max_attempts_exceeded

  def _finished(self):
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
      for foo in itertbx.islice(flipping,
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
        observable = skewness_evolution()
        for n, flipping in enumerate(
          itertbx.islice(self.flipping_iterator,
                         0, self.max_solving_iterations)):
          self.iteration_index = n
          if n % self.yield_solving_interval == 0:
            yield self.solving
          observable.append(flipping.rho_map.skewness())
          if flipping.rho_map.skewness() < 3: continue
          if observable.had_phase_transition():
            self.attempts.append(n)
            yield self.polishing
            break
        else:
          if i_attempt != self.max_attempts_to_get_phase_transition:
            self.flipping_iterator.restart()
            yield self.guessing_delta
      self.max_attempts_exceeded = True
      yield self.finished

  def _polishing(self):
    while 1:
      self.flipping_iterator.__class__ = low_density_elimination_iterator
      self.flipping_iterator.constant_rho_c = self.flipping_iterator.delta
      self.flipping_iterator.restart(
        f_000=0,
        f_calc=self.flipping_iterator.f_calc)
      for i in xrange(self.polishing_iterations):
        self.flipping_iterator.next()
      yield self.evaluating

  def _evaluating(self):
    while 1:
      attempts = 0
      while attempts < self.max_attempts_to_get_sharp_correlation_map:
        attempts += 1
        self.f_calc_solutions = []
        for f_calc, shift, cc_peak_height\
                           in self.flipping_iterator.f_calc_symmetrisations():
          if not self.f_calc_solutions and cc_peak_height < 0.9:
            yield self.guessing_delta
            break
          self.f_calc_solutions.append((f_calc, shift, cc_peak_height))
        else:
          yield self.finished
      selt.max_attempts_exceeded = True

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
        print "Polishing with rho_c=%.4f" % flipping.rho_c()
    elif solving.state is solving.finished:
      break

    previous_state = solving.state



class skewness_evolution(object):

  smoothing_coefficient = 0.25
  increasing = True
  noise_level_before = 0.3
  noise_level_after = 0.2

  def __init__(self, **kwds):
    adopt_optional_init_args(self, kwds)
    self.values = flex.double()
    self.differences = flex.double()

  def append(self, x):
    if len(self.values) > 1:
      y0 = self.values[-1]
      a = self.smoothing_coefficient
      y1 = y0 + a*(x - y0)
      self.values.append(y1)
      delta = y1 - y0
      if not self.increasing: delta = -delta
      self.differences.append(delta)
    else:
      self.values.append(x)

  def had_phase_transition(self):
    if len(self.differences) < 5: return False
    i_max = flex.max_index(self.differences)
    noise_before = (self.differences
                    < self.noise_level_before*self.differences[i_max])
    before = flex.last_index(noise_before[:i_max], True)
    if before is None: before = -1
    before += 1
    if i_max - before < 4: return False
    negative_after = self.differences < 0
    after = flex.first_index(negative_after[i_max:], True)
    if after is None: return False
    after += i_max
    if after - before < 10: return False
    if len(self.values) - after < 10: return False
    tail_stats = scitbx.math.basic_statistics(self.differences[-5:])
    if (tail_stats.max_absolute
        > self.noise_level_after*self.differences[i_max]): return False
    return True
