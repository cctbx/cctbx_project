""" Tools to solve non-linear L.S. problems formulated with normal-equations.
"""

import libtbx
from scitbx.array_family import flex
import itertools

class iterations(object):

  track_shifts = False
  track_gradients = False
  track_all = False

  def __init__(self, normal_eqns, **kwds):
    """ normal_eqns, thereafter denoted as e, shall feature:
          - e.build_up():
              build the normal equations
          - e.objective:
              the value of the L.S. objective
          - e.gradient:
              the gradient of the objective wrt the independent parameters
          - e.gradient_norm_reference:
              a scalar to compare the norm of the gradient to

        An instance of this class accumulates the following information:
        L.S. objective and gradient, relative gradient norm
        (i.e. gradient norm over the gradient norm reference),
        and L.S. scale factor if the normal equations feature one.
    """
    self.normal_eqns = normal_eqns
    libtbx.adopt_optional_init_args(self, kwds)
    self.objectives = flex.double()
    self.gradients = []
    self.relative_gradient_norms = flex.double()
    self.shifts = []
    if hasattr(normal_eqns, "scale_factor"):
      self.scale_factors = flex.double()
    else:
      self.scale_factors = None

  def accumulate_normal_equation_properties(self):
    self.objectives.append(self.normal_eqns.objective)
    g = self.normal_eqns.gradient
    self.relative_gradient_norms.append(
      g.norm()/self.normal_eqns.gradient_norm_reference)
    if self.scale_factors is not None:
      self.scale_factors.append(self.normal_eqns.scale_factor)
    if self.track_all or self.track_gradients:
      self.gradients.append(self.normal_eqns.gradient)

  def accumulate_solution_properties(self):
    if self.track_all or self.track_shifts:
      self.shifts.append(self.normal_eqns.shifts)

  def accumulate(self):
    self.accumulate_normal_equation_properties()
    self.accumulate_solution_properties()

  def do(self,
         relative_gradient_norm_threshold=None,
         n_iterations=None):
    assert (n_iterations is not None
            or relative_gradient_norm_threshold is not None)
    if n_iterations is None:
      counter = itertools.count()
    else:
      counter = xrange(n_iterations)
    for i in counter:
      self.next()
      if relative_gradient_norm_threshold is None: continue
      if self.relative_gradient_norms[-1] < relative_gradient_norm_threshold:
        break

  class n_iterations(libtbx.property):
    def fget(self): return len(self.objectives)

class naive_iterations(iterations):

  def next(self):
    self.normal_eqns.build_up()
    self.accumulate_normal_equation_properties()
    self.normal_eqns.solve_and_apply_shifts()
    self.accumulate_solution_properties()
    return self
