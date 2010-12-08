from __future__ import division

""" Tools to solve non-linear L.S. problems formulated with normal-equations.
"""

import libtbx
from scitbx.array_family import flex


class journaled_normal_eqns(object):
  """ A decorator that keeps the history of the objective, gradient,
  shifts, etc of an underlying normal equations object. An instance of this
  class is a drop-in replacement of that underlying object, with the
  journaling just mentionned automatically happening behind the scene.
  """

  def __init__(self, normal_eqns, journal, track_gradient, track_shifts):
    """ Decorate the given normal equations. The history will be accumulated
    in relevant attributes of journal. The flags track_xxx specify whether
    to journal the gradient and/or the shifts, a potentially memory-hungry
    operation.
    """
    self.actual = normal_eqns
    self.journal = journal
    self.journal.objective_history = flex.double()
    if track_gradient:
      self.journal.gradient_history = []
    else:
      self.journal.gradient_history = None
    self.journal.gradient_norm_history = flex.double()
    if track_shifts:
      self.journal.shifts_history = []
    else:
      self.journal.shifts_history = None
    self.journal.shifts_norm_history = flex.double()
    self.journal.parameter_vector_norm_history = flex.double()
    if hasattr(normal_eqns, "scale_factor"):
      self.journal.scale_factor_history = flex.double()
    else:
      self.journal.scale_factor_history = None

  def __getattr__(self, name):
    return getattr(self.actual, name)

  def build_up(self):
    self.actual.build_up()
    self.journal.parameter_vector_norm_history.append(
      self.actual.parameter_vector_norm)
    self.journal.objective_history.append(self.actual.objective)
    self.journal.gradient_norm_history.append(self.actual.gradient.norm_inf())
    if self.journal.gradient_history is not None:
      self.journal.gradient_history.append(self.actual.gradient)
    if self.journal.scale_factor_history is not None:
      self.journal.scale_factor_history.append(self.actual.scale_factor)

  def solve(self):
    self.actual.solve()
    self.journal.shifts_norm_history.append(self.actual.shifts.norm())
    if self.journal.shifts_history is not None:
      self.journal.shifts_history.append(self.actual.shifts.deep_copy())

  def apply_shifts(self):
    self.actual.apply_shifts()

  def solve_and_apply_shifts(self):
    self.solve()
    self.apply_shifts()


class iterations(object):
  """ Iterations to solve a L.S. minimisation problem.

  It is assumed that the objective function is properly scaled to be
  between 0 and 1 (or approximately so), so that this class can
  meaningfully thresholds the norm of the gradients to watch for
  convergence.

  Use classic stopping criteria: c.f. e.g.
  Methods for non-linear least-squares problems,
  K. Madsen, H.B. Nielsen, O. Tingleff,
  http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
  """

  track_shifts = False
  track_gradients = False
  track_all = False
  n_max_iterations = 100
  gradient_threshold = None
  shift_threshold = None

  def __init__(self, normal_eqns, **kwds):
    """
    """
    libtbx.adopt_optional_init_args(self, kwds)
    if self.track_all: self.track_shifts = self.track_gradients = True
    self.normal_eqns = journaled_normal_eqns(normal_eqns, self,
                                             self.track_gradients,
                                             self.track_shifts)
    self.do()

  def has_gradient_converged_to_zero(self):
    eps_1 = self.gradient_threshold
    return eps_1 is not None and self.gradient_norm_history[-1] <= eps_1

  def had_too_small_a_step(self):
    eps_2 = self.shift_threshold
    if eps_2 is None: return False
    x = self.parameter_vector_norm_history[-1]
    h = self.shifts_norm_history[-1]
    return h <= eps_2*(x + eps_2)

  def do(self):
    raise NotImplementedError


class naive_iterations(iterations):

  def do(self):
    self.n_iterations = 0
    while self.n_iterations <= self.n_max_iterations:
      self.normal_eqns.build_up()
      if self.has_gradient_converged_to_zero(): break
      self.normal_eqns.solve()
      if self.had_too_small_a_step(): break
      self.normal_eqns.apply_shifts()
      self.n_iterations += 1

  def __str__(self):
    return "pure Gauss-Newton"
