""" Tools to solve non-linear L.S. problems formulated with normal-equations.
"""

import libtbx
from scitbx.array_family import flex
import itertools

class iterations(object):

  track_shifts = False
  track_gradients = False
  track_all = False
  n_max_iterations = 100

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
    self.gradient_norms.append(g.norm_1())
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
