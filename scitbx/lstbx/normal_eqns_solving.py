from __future__ import division

""" Tools to solve non-linear L.S. problems formulated with normal-equations.
"""

import libtbx
from scitbx.array_family import flex


class journaled_non_linear_ls(object):
  """ A decorator that keeps the history of the objective, gradient,
  step, etc of an underlying normal equations object. An instance of this
  class is a drop-in replacement of that underlying object, with the
  journaling just mentioned automatically happening behind the scene.
  """

  def __init__(self, non_linear_ls, journal, track_gradient, track_step):
    """ Decorate the given normal equations. The history will be accumulated
    in relevant attributes of journal. The flags track_xxx specify whether
    to journal the gradient and/or the step, a potentially memory-hungry
    operation.
    """
    self.actual = non_linear_ls
    self.journal = journal
    self.journal.objective_history = flex.double()
    if track_gradient:
      self.journal.gradient_history = []
    else:
      self.journal.gradient_history = None
    self.journal.gradient_norm_history = flex.double()
    if track_step:
      self.journal.step_history = []
    else:
      self.journal.step_history = None
    self.journal.step_norm_history = flex.double()
    self.journal.parameter_vector_norm_history = flex.double()
    if hasattr(non_linear_ls, "scale_factor"):
      self.journal.scale_factor_history = flex.double()
    else:
      self.journal.scale_factor_history = None

  def __getattr__(self, name):
    return getattr(self.actual, name)

  def build_up(self, objective_only=False):
    self.actual.build_up(objective_only)
    if objective_only: return
    self.journal.parameter_vector_norm_history.append(
      self.actual.parameter_vector_norm())
    self.journal.objective_history.append(self.actual.objective())
    self.journal.gradient_norm_history.append(
      self.actual.opposite_of_gradient().norm_inf())
    if self.journal.gradient_history is not None:
      self.journal.gradient_history.append(-self.actual.opposite_of_gradient())
    if self.journal.scale_factor_history is not None:
      self.journal.scale_factor_history.append(self.actual.scale_factor())

  def solve(self):
    self.actual.solve()
    self.journal.step_norm_history.append(self.actual.step().norm())
    if self.journal.step_history is not None:
      self.journal.step_history.append(self.actual.step().deep_copy())

  def step_forward(self):
    self.actual.step_forward()

  def solve_and_step_forward(self):
    self.solve()
    self.step_forward()


class iterations(object):
  """ Iterations to solve a non-linear L.S. minimisation problem.

  It is assumed that the objective function is properly scaled to be
  between 0 and 1 (or approximately so), so that this class can
  meaningfully thresholds the norm of the gradients to watch for
  convergence.

  The interface expected from the normal equations object passed to __init__
  is that of lstbx.non_linear_normal_equations_mixin

  Use classic stopping criteria: c.f. e.g.
  Methods for non-linear least-squares problems,
  K. Madsen, H.B. Nielsen, O. Tingleff,
  http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf

  The do_damping and do_scale_shifts function implements the shelxl damping as
  described here:
  http://shelx.uni-ac.gwdg.de/SHELX/shelx.pdf
  do_scale_shifts aslo compares the maximum shift/esd with
  convergence_as_shift_over_esd and returns boolean to indicate if
  the refinement converged according to this criterion
  These functions do nothing unless called by a derived class
  """

  track_step = False
  track_gradient = False
  track_all = False
  n_max_iterations = 100
  gradient_threshold = None
  step_threshold = None
  damping_value = 0.0007
  max_shift_over_esd = 15
  convergence_as_shift_over_esd = 1e-5

  def __init__(self, non_linear_ls, **kwds):
    """
    """
    libtbx.adopt_optional_init_args(self, kwds)
    if self.track_all: self.track_step = self.track_gradient = True
    self.non_linear_ls = journaled_non_linear_ls(non_linear_ls, self,
                                                 self.track_gradient,
                                                 self.track_step)
    self.do()

  def has_gradient_converged_to_zero(self):
    eps_1 = self.gradient_threshold
    g = self.gradient_norm_history[-1]
    return eps_1 is not None and g <= eps_1

  def had_too_small_a_step(self):
    eps_2 = self.step_threshold
    if eps_2 is None: return False
    x = self.parameter_vector_norm_history[-1]
    h = self.step_norm_history[-1]
    return h <= eps_2*(x + eps_2)

  def do_damping(self, value):
    a = self.non_linear_ls.normal_matrix_packed_u()
    a.matrix_packed_u_diagonal_add_in_place(value*a.matrix_packed_u_diagonal())

  def do_scale_shifts(self, max_shift_over_esd):
    x = self.non_linear_ls.step()
    esd = self.non_linear_ls.covariance_matrix().matrix_packed_u_diagonal()
    x_over_esd = flex.abs(x/flex.sqrt(esd))
    max_val = flex.max(x_over_esd)
    if max_val < self.convergence_as_shift_over_esd:
      return True
    if max_val > max_shift_over_esd:
      shift_scale = max_shift_over_esd/max_val
      x *= shift_scale
    return False

  def do(self):
    raise NotImplementedError


class naive_iterations(iterations):

  def do(self):
    self.n_iterations = 0
    while self.n_iterations < self.n_max_iterations:
      self.non_linear_ls.build_up()
      if self.has_gradient_converged_to_zero(): break
      self.non_linear_ls.solve()
      if self.had_too_small_a_step(): break
      self.non_linear_ls.step_forward()
      self.n_iterations += 1

  def __str__(self):
    return "pure Gauss-Newton"


class naive_iterations_with_damping(iterations):

  def do(self):
    self.n_iterations = 0
    do_last = False
    while self.n_iterations < self.n_max_iterations:
      self.non_linear_ls.build_up()
      if self.has_gradient_converged_to_zero():
        do_last = True
      if not do_last and self.n_iterations+1 < self.n_max_iterations:
        self.do_damping(self.damping_value)
      self.non_linear_ls.solve()
      step_too_small = self.had_too_small_a_step()
      self.non_linear_ls.step_forward()
      self.n_iterations += 1
      if do_last: break
      if step_too_small:
        do_last = True

  def __str__(self):
    return "pure Gauss-Newton with damping"


class naive_iterations_with_damping_and_shift_limit(iterations):

  def do(self):
    self.n_iterations = 0
    do_last = False
    while self.n_iterations < self.n_max_iterations:
      self.non_linear_ls.build_up()
      if self.has_gradient_converged_to_zero():
        do_last = True
      if not do_last and self.n_iterations+1 < self.n_max_iterations:
        self.do_damping(self.damping_value)
      self.non_linear_ls.solve()
      step_too_small = self.had_too_small_a_step()
      if not do_last and not step_too_small:
        step_too_small = self.do_scale_shifts(self.max_shift_over_esd)
      self.non_linear_ls.step_forward()
      self.n_iterations += 1
      if do_last: break
      if step_too_small:
        do_last = True

  def __str__(self):
    return "pure Gauss-Newton with damping and shift scaling"


class levenberg_marquardt_iterations(iterations):

  tau = 1e-3

  class mu(libtbx.property):
    def fget(self):
      return self._mu
    def fset(self, value):
      self.mu_history.append(value)
      self._mu = value

  def do(self):
    self.mu_history = flex.double()
    self.n_iterations = 0
    nu = 2
    self.non_linear_ls.build_up()
    if self.has_gradient_converged_to_zero(): return
    a = self.non_linear_ls.normal_matrix_packed_u()
    self.mu = self.tau*flex.max(a.matrix_packed_u_diagonal())
    while self.n_iterations < self.n_max_iterations:
      a.matrix_packed_u_diagonal_add_in_place(self.mu)
      objective = self.non_linear_ls.objective()
      g = -self.non_linear_ls.opposite_of_gradient()
      self.non_linear_ls.solve()
      if self.had_too_small_a_step(): break
      self.n_iterations += 1
      h = self.non_linear_ls.step()
      expected_decrease = 0.5*h.dot(self.mu*h - g)
      self.non_linear_ls.step_forward()
      self.non_linear_ls.build_up(objective_only=True)
      objective_new = self.non_linear_ls.objective()
      actual_decrease = objective - objective_new
      rho = actual_decrease/expected_decrease
      if rho > 0:
        if self.has_gradient_converged_to_zero(): break
        self.mu *= max(1/3, 1 - (2*rho - 1)**3)
        nu = 2
      else:
        self.non_linear_ls.step_backward()
        self.mu *= nu
        nu *= 2
      self.non_linear_ls.build_up()

  def __str__(self):
    return "Levenberg-Marquardt"
