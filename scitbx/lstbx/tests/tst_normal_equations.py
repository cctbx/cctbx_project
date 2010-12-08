from scitbx.array_family import flex
from scitbx import sparse
from scitbx.lstbx import normal_eqns, normal_eqns_solving
from libtbx.test_utils import approx_equal, Exception_expected
import libtbx
import itertools

def exercise_linear_normal_equations():
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


class polynomial_fit(normal_eqns.non_linear_ls_with_separable_scale_factor,
                     normal_eqns.non_linear_ls_mixin):

  noise = 1e-5
  n_data = 10
  x_0 = flex.double((0.5, 0.3, 0.2))

  def __init__(self, **kwds):
    super(polynomial_fit, self).__init__(n_parameters=3, normalised=False)
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

def exercise_non_linear_ls_with_separable_scale_factor():
  test = polynomial_fit()
  test.build_up()
  # Reference values computed in tst_normal_equations.nb
  eps = 5e-14
  assert approx_equal(test.optimal_scale_factor(), 0.6148971786833856, eps)
  assert approx_equal(test.objective(), 0.039642707534326034, eps)


  assert not test.step_equations().solved
  try:
    test.step_equations().cholesky_factor_packed_u()
    raise Exception_expected
  except RuntimeError:
    pass
  try:
    test.step_equations().solution()
    raise Exception_expected
  except RuntimeError:
    pass

  assert approx_equal(
    list(test.step_equations().normal_matrix_packed_u()),
    [ 0.371944193675858, 0.39066546997866547  , 0.10797294655500618,
                           0.41859250354804045, 0.08077629438075473,
                                                0.19767268057900367 ],
    eps)
  assert approx_equal(
    list(test.step_equations().right_hand_side()),
    [ 0.12149917297914861, 0.13803759252793774, -0.025190641142579157 ],
    eps)

  test.step_equations().solve()

  assert test.step_equations().solved
  try:
    test.step_equations().normal_matrix_packed_u()
    raise Exception_expected
  except RuntimeError:
    pass
  try:
    test.step_equations().right_hand_side()
    raise Exception_expected
  except RuntimeError:
    pass

  assert approx_equal(
    list(test.step_equations().cholesky_factor_packed_u()),
    [ 0.6098722765266986, 0.6405693208478925   ,  0.1770418999366983 ,
                            0.09090351333425013, -0.3589664912436558 ,
                                                  0.19357661121640218 ],
    eps)
  assert approx_equal(
    list(test.step_equations().solution()),
    [ 1.2878697604109028, -0.7727798877778043, -0.5151113342942297 ],
    eps=1e-12)


class exponential_fit(
  normal_eqns.non_linear_ls,
  normal_eqns.non_linear_ls_mixin):

  """ Model M(x, t) = x3 e^{x1 t} + x4 e^{x2 t} """

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
  exercise_linear_normal_equations()
  exercise_non_linear_ls_with_separable_scale_factor()
  print 'OK'

if __name__ == '__main__':
  run()
