""" A collection of L.S. test problems """

from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns
import libtbx

class polynomial_fit(normal_eqns.non_linear_ls_with_separable_scale_factor,
                     normal_eqns.non_linear_ls_mixin):
  """ Fit yc(t) = -t^3 + a t^2 + b t + c to synthetic data produced by
      yo(t) = 2 t^2 (1-t) + noise, with an overall scale factor on yc(t)
      whose optimal value is therefore 2.
  """

  noise = 1e-5
  n_data = 10
  x_0 = flex.double((0.5, 0.3, 0.2))

  def __init__(self, normalised, **kwds):
    super(polynomial_fit, self).__init__(n_parameters=3,
                                         normalised=normalised)
    libtbx.adopt_optional_init_args(self, kwds)
    self.t = t = flex.double_range(self.n_data)/self.n_data
    noise = self.noise*flex.double([ (-1)**i for i in xrange(self.n_data) ])
    self.yo = 2.*t**2.*(1-t) + noise
    self.one = flex.double(self.n_data, 1)
    self.t2 = t**2
    self.t3 = t**3
    self.restart()

  def restart(self):
    self.x = self.x_0.deep_copy()
    self.old_x = None

  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None

  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    a, b, c = self.x
    one, t, t2, t3 = self.one, self.t, self.t2, self.t3
    yc = -t3 + a*t2 + b*t + c
    grad_yc = (t2, t, one)
    jacobian_yc = flex.double(flex.grid(self.n_data, self.n_parameters))
    for j, der_r in enumerate(grad_yc):
      jacobian_yc.matrix_paste_column_in_place(der_r, j)
    self.add_equations(yc, jacobian_yc, self.yo, weights=None)
    self.finalise()


class polynomial_fit_with_penalty(polynomial_fit):

  def build_up(self, objective_only=False):
    super(polynomial_fit_with_penalty, self).build_up(objective_only)
    a, b, c = self.x
    reduced = self.reduced_problem()
    reduced.add_equation(residual=(a - b + c - 2),
                         grad_residual=flex.double((1, -1, 1)),
                         weight=1)

class exponential_fit(
  normal_eqns.non_linear_ls,
  normal_eqns.non_linear_ls_mixin):

  """ Model M(x, t) = x3 e^{x1 t} + x4 e^{x2 t}

      Problem 18 from

      UCTP Test Problems for Unconstrained Optimization
      Hans Bruun Nielsen
      TECHNICAL REPORT IMM-REP-2000-17
  """

  n_data = 45
  arg_min = flex.double((-4, -5, 4, -4))
  x_0     = flex.double((-1, -2, 1, -1))

  def __init__(self):
    super(exponential_fit, self).__init__(n_parameters=4)
    self.t = 0.02*flex.double_range(1, self.n_data + 1)
    self.y = flex.double((
      0.090542, 0.124569, 0.179367, 0.195654, 0.269707, 0.286027, 0.289892,
      0.317475, 0.308191, 0.336995, 0.348371, 0.321337, 0.299423, 0.338972,
      0.304763, 0.288903, 0.300820, 0.303974, 0.283987, 0.262078, 0.281593,
      0.267531, 0.218926, 0.225572, 0.200594, 0.197375, 0.182440, 0.183892,
      0.152285, 0.174028, 0.150874, 0.126220, 0.126266, 0.106384, 0.118923,
      0.091868, 0.128926, 0.119273, 0.115997, 0.105831, 0.075261, 0.068387,
      0.090823, 0.085205, 0.067203
      ))
    assert len(self.y) == len(self.t)
    self.restart()

  def restart(self):
    self.x = self.x_0.deep_copy()
    self.old_x = None

  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    x1, x2, x3, x4 = self.x
    exp_x1_t = flex.exp(x1*self.t)
    exp_x2_t = flex.exp(x2*self.t)
    residuals = x3*exp_x1_t + x4*exp_x2_t
    residuals -= self.y

    self.reset()
    if objective_only:
      self.add_residuals(residuals, weights=None)
    else:
      grad_r = (self.t*x3*exp_x1_t,
                self.t*x4*exp_x2_t,
                exp_x1_t,
                exp_x2_t)
      jacobian = flex.double(flex.grid(self.n_data, self.n_parameters))
      for j, der_r in enumerate(grad_r):
        jacobian.matrix_paste_column_in_place(der_r, j)
      self.add_equations(residuals, jacobian, weights=None)

  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None
