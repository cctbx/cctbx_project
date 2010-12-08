from scitbx.array_family import flex
from scitbx import sparse
from scitbx.lstbx import normal_eqns, normal_eqns_solving
from libtbx.test_utils import approx_equal
import itertools

def exercise_basic_normal_equations():
  py_eqs = [ ( 1, (-1,  0,  0),  1),
             ( 2, ( 2, -1,  0),  3),
             (-1, ( 0,  2,  1),  2),
             (-2, ( 0,  1,  0), -2),
             ]

  eqs_0 = normal_eqns.linear_ls(3)
  for b, a, w in py_eqs:
    eqs_0.add_equation(right_hand_side=b,
                       design_matrix_row=flex.double(a),
                       weight=w)

  eqs_1 = normal_eqns.linear_ls(3)
  b = flex.double()
  w = flex.double()
  a = sparse.matrix(len(py_eqs), 3)
  for i, (b_, a_, w_) in enumerate(py_eqs):
    b.append(b_)
    w.append(w_)
    for j in xrange(3):
      if a_[j]: a[i, j] = a_[j]
  eqs_1.add_equations(right_hand_side=b, design_matrix=a, weights=w)

  assert approx_equal(
    eqs_0.normal_matrix_packed_u(), eqs_1.normal_matrix_packed_u(), eps=1e-15)
  assert approx_equal(
    eqs_0.right_hand_side(), eqs_1.right_hand_side(), eps=1e-15)
  assert approx_equal(
    list(eqs_0.normal_matrix_packed_u()), [ 13, -6, 0, 9, 4, 2 ], eps=1e-15)
  assert approx_equal(
    list(eqs_0.right_hand_side()), [ 11, -6, -2 ], eps=1e-15)

def exercise_normal_equations_separating_scale_factor():
  eqs = lstbx.non_linear_ls_with_separable_scale_factor(3)
  eqs.add_equation(y_calc=1.1,
                   grad_y_calc=flex.double((1, 2, 3)),
                   y_obs=1,
                   weight=1)
  eqs.add_equation(y_calc=2.2,
                   grad_y_calc=flex.double((2, 3, 1)),
                   y_obs=2,
                   weight=1)
  eqs.add_equation(y_calc=3.3,
                   grad_y_calc=flex.double((3, 1, 2)),
                   y_obs=3,
                   weight=1)
  eqs.add_equation(y_calc=4.4,
                   grad_y_calc=flex.double((1, 3, 2)),
                   y_obs=4,
                   weight=1)
  eqs.finalise()
  a, b = eqs.step_equations()
  assert a.size() == 6
  assert b.size() == 3
  eqs.step_equations().solve()
  assert eqs.step_equations().solved



class exponential_fit(
  normal_eqns.non_linear_ls,
  normal_eqns.non_linear_ls_mixin):

  """ Model M(x, t) = x3 e^{x1 t} + x4 e^{x2 t} """

  arg_min = flex.double((-4, -5, 4, -4))
  x_0     = flex.double((-1, -2, 1, -1))

  def __init__(self):
    super(exponential_fit, self).__init__(n_parameters=4)
    self.t = 0.02*flex.double_range(1, 46)
    self.y = flex.double((
      0.090542, 0.124569, 0.179367, 0.195654, 0.269707, 0.286027, 0.289892,
      0.317475, 0.308191, 0.336995, 0.348371, 0.321337, 0.299423, 0.338972,
      0.304763, 0.288903, 0.300820, 0.303974, 0.283987, 0.262078, 0.281593,
      0.267531, 0.218926, 0.225572, 0.200594, 0.197375, 0.182440, 0.183892,
      0.152285, 0.174028, 0.150874, 0.126220, 0.126266, 0.106384, 0.118923,
      0.091868, 0.128926, 0.119273, 0.115997, 0.105831, 0.075261, 0.068387,
      0.090823, 0.085205, 0.067203
      ))
    self.restart()

  def restart(self):
    self.x = self.x_0.deep_copy()
    self.old_x = None
    self.old_objective = None

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
      jacobian = flex.double(flex.grid(len(self.t), len(self.x)))
      for j, der_r in enumerate(grad_r):
        jacobian.matrix_paste_column_in_place(der_r, j)
      self.add_equations(residuals, jacobian, weights=None)

  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None


def exercise_levenberg_marquardt(non_linear_ls, plot=False):
  non_linear_ls.restart()
  iterations = normal_eqns_solving.levenberg_marquardt_iterations(
    non_linear_ls,
    track_all=True,
    gradient_threshold=1e-8,
    step_threshold=1e-8,
    tau=1e-4,
    n_max_iterations=200)
  assert approx_equal(non_linear_ls.x, non_linear_ls.arg_min, eps=5e-4)
  print "L-M: %i iterations" % iterations.n_iterations
  if plot:
    f = open('plot.nb', 'w')
    print >>f, "g=%s;" % iterations.gradient_norm_history.mathematica_form()
    print >>f, "\[Mu]=%s;" % iterations.mu_history.mathematica_form()
    print >>f, "ListLogPlot[{g,\[Mu]},Joined->True]"
    f.close()

def run():
  import sys
  plot = '--plot' in sys.argv[1:]
  t = exponential_fit()
  exercise_levenberg_marquardt(t, plot)
  exercise_basic_normal_equations()
  exercise_normal_equations_separating_scale_factor()
  print 'OK'

if __name__ == '__main__':
  run()
