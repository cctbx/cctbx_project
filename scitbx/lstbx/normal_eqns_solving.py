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
    if hasattr(normal_eqns, "scale_factor"):
      self.journal.scale_factor_history = flex.double()
    else:
      self.journal.scale_factor_history = None

  def __getattr__(self, name):
    return getattr(self.actual, name)

  def build_up(self):
    self.actual.build_up()
    self.journal.objective_history.append(self.actual.objective)
    self.journal.gradient_norm_history.append(self.actual.gradient.norm_inf())
    if self.journal.gradient_history is not None:
      self.journal.gradient_history.append(self.actual.gradient)
    if self.journal.scale_factor_history is not None:
      self.journal.scale_factor_history.append(self.actual.scale_factor)

  def solve(self):
    self.actual.solve()
    if self.journal.shifts_history is not None:
      self.journal.shifts_history.append(self.actual.shifts.deep_copy())

  def apply_shifts(self):
    self.actual.apply_shifts()


class iterations(object):

  track_shifts = False
  track_gradients = False
  track_all = False
  n_max_iterations = 100
  norm_of_gradient = flex.double.norm_inf

  def __init__(self, normal_eqns, **kwds):
    """ normal_eqns, thereafter denoted as e, shall feature:
          - e.build_up():
              build the normal equations
          - e.objective:
              the value of the L.S. objective
          - e.gradient:
              the gradient of the objective wrt the independent parameters

        It is assumed that the objective is properly scaled to be
        between 0 and 1 (or approximately so), so that this class can
        meaningfully thresholds the norm of the gradients to watch for
        convergence.

        An instance of this class accumulates the following information:
        L.S. objective and gradient, and L.S. scale factor
        if the normal equations feature one.
    """
    self.normal_eqns = normal_eqns
    libtbx.adopt_optional_init_args(self, kwds)
    self.objectives = flex.double()
    self.gradients = []
    self.gradient_norms = []
    self.shifts = []
    if hasattr(normal_eqns, "scale_factor"):
      self.scale_factors = flex.double()
    else:
      self.scale_factors = None

  class n_iterations(libtbx.property):
    def fget(self): return len(self.objectives)

  def accumulate_normal_equation_properties(self):
    self.objectives.append(self.normal_eqns.objective)
    g = self.normal_eqns.gradient
    self.gradient_norms.append(self.norm_of_gradient(g))
    if self.scale_factors is not None:
      self.scale_factors.append(self.normal_eqns.scale_factor)
    if self.track_all or self.track_gradients:
      self.gradients.append(g)

  def accumulate_solution_properties(self):
    if self.track_all or self.track_shifts:
      self.shifts.append(self.normal_eqns.shifts.deep_copy())

  def accumulate(self):
    self.accumulate_normal_equation_properties()
    self.accumulate_solution_properties()

  def do(self,
         gradient_threshold=None,
         shift_threshold=None,
         n_iterations = None):
    """ Iterate, using classic stopping criteria (c.f. e.g.
        Methods for non-linear least-squares problems,
        K. Madsen, H.B. Nielsen, O. Tingleff,
        http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
    """
    eps_1 = gradient_threshold
    eps_2 = shift_threshold
    for i in xrange(n_iterations or self.n_max_iterations):
      self.next()
      h = self.normal_eqns.shifts
      if (eps_1 is not None and self.gradient_norms[-1] <= eps_1):
        break
      if eps_2 is not None:
        lxl = self.normal_eqns.parameter_vector_norm
        if h.norm() <= eps_2*(lxl + eps_2):
          break


class naive_iterations(iterations):

  def next(self):
    self.normal_eqns.build_up()
    self.accumulate_normal_equation_properties()
    self.normal_eqns.solve_and_apply_shifts()
    self.accumulate_solution_properties()
    return self
