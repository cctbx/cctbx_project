from __future__ import absolute_import, division, print_function
from scitbx.examples import immoptibox_ports
from scitbx.math import gaussian
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, eps_eq
from libtbx.utils import format_cpu_times
from six.moves import range
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle
from six.moves import cStringIO as StringIO
import math
import sys

def finite_gradient_dx_at_x(gaussian, x, eps=1.e-5):
  if (x == 0): return 0
  assert x >= eps
  tm = gaussian.at_x(x-eps)
  tp = gaussian.at_x(x+eps)
  return (tp-tm)/(2*eps)

def exercise_gradient_dx(gaussian, x_max=1., n_points=50):
  for i in range(n_points+1):
    x = x_max * i / n_points
    grad_finite = finite_gradient_dx_at_x(gaussian, x)
    grad_analytical = gaussian.gradient_dx_at_x(x)
    assert eps_eq(grad_finite, grad_analytical)

def exercise_integral_dx(gaussian, x_max=1., n_points=1000):
  numerical_integral = 0
  x_step = x_max / n_points
  for i in range(n_points+1):
    x = x_max * i / n_points
    new_value = gaussian.at_x(x)
    if (i):
      numerical_integral += (prev_value + new_value) * .5
    prev_value = new_value
    analytical_integral = gaussian.integral_dx_at_x(x, 1.e-3)
    assert eps_eq(analytical_integral, gaussian.integral_dx_at_x(x))
    assert eps_eq(numerical_integral*x_step, analytical_integral, eps=1.e-5)

def term_finite_gradient_d_ab_at_x(term, x, eps=1.e-5):
  tm = gaussian.term(term.a-eps,term.b).at_x(x)
  tp = gaussian.term(term.a+eps,term.b).at_x(x)
  gr_a = (tp-tm)/(2*eps)
  tm = gaussian.term(term.a,term.b-eps).at_x(x)
  tp = gaussian.term(term.a,term.b+eps).at_x(x)
  gr_b = (tp-tm)/(2*eps)
  return gaussian.term(gr_a, gr_b)

def exercise_term_gradients_d_ab(term, x_max=1., n_points=50):
  for i in range(n_points+1):
    x = x_max * i / n_points
    grad_finite = term_finite_gradient_d_ab_at_x(term, x)
    grad_analytical = term.gradients_d_ab_at_x_sq(x*x)
    assert eps_eq(grad_finite.a, grad_analytical.a)
    assert eps_eq(grad_finite.b, grad_analytical.b)

def exercise_term():
  t = gaussian.term(2,3)
  assert approx_equal(t.a, 2)
  assert approx_equal(t.b, 3)
  assert approx_equal(t.at_x_sq(4), 2*math.exp(-3*4))
  assert approx_equal(t.at_x(2), 2*math.exp(-3*4))
  eps = 1.e-5
  for ix in (range(10)):
    x = ix/10.
    assert eps_eq((t.at_x(x+eps)-t.at_x(x-eps))/(2*eps), t.gradient_dx_at_x(x))
  for f in [1,-1]:
    for t in [gaussian.term(f*2,3),
              gaussian.term(f*3,0),
              gaussian.term(f*4,1.e-4),
              gaussian.term(f*5,-1)]:
      exercise_gradient_dx(t)
      exercise_integral_dx(t)
      exercise_term_gradients_d_ab(t)

def exercise_sum():
  g = gaussian.sum(0)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  assert g.n_parameters() == 1
  assert approx_equal(g.parameters(), [0])
  g = gaussian.sum(0, True)
  assert g.use_c()
  g = gaussian.sum(0, False)
  assert not g.use_c()
  g = gaussian.sum(1)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), 1)
  assert g.use_c()
  assert g.n_parameters() == 1
  assert approx_equal(g.parameters(), [1])
  g = gaussian.sum((), ())
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert g.c() == 0
  assert not g.use_c()
  assert g.n_parameters() == 0
  assert g.parameters().size() == 0
  g = gaussian.sum((), (), -2)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), -2)
  g = gaussian.sum(flex.double((1,2,3,4)))
  assert approx_equal(g.array_of_a(), (1,3))
  assert approx_equal(g.array_of_b(), (2,4))
  assert approx_equal(g.c(), 0)
  assert not g.use_c()
  assert approx_equal(g.parameters(), [1,2,3,4])
  g = gaussian.sum(flex.double((1,2,3,4)), 0, True)
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  g = gaussian.sum(flex.double((1,2,3,4)), 5)
  assert approx_equal(g.c(), 5)
  assert g.use_c()
  assert approx_equal(g.parameters(), [1,2,3,4,5])
  g = gaussian.sum(flex.double((1,2,3,4,5)))
  assert approx_equal(g.c(), 5)
  assert g.use_c()
  assert approx_equal(g.parameters(), [1,2,3,4,5])
  g = gaussian.sum((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  assert g.n_terms() == 5
  assert approx_equal(g.array_of_a(),(1,-2,3,-4,5))
  assert approx_equal(g.array_of_b(),(-.1,.2,-.3,.4,-.5))
  assert approx_equal(g.c(), 6)
  assert approx_equal(g.at_x_sq(3/4.), 13.4251206)
  assert approx_equal(g.at_x_sq(flex.double([2/4.,3/4.])),
    [11.8723031, 13.4251206])
  assert approx_equal(g.at_x(math.sqrt(3/4.)), 13.4251206)
  assert approx_equal(g.at_x(flex.sqrt(flex.double([2/4.,3/4.]))),
    [11.8723031, 13.4251206])
  s = pickle.dumps(g)
  l = pickle.loads(s)
  assert l.n_terms() == g.n_terms()
  assert approx_equal(l.array_of_a(), g.array_of_a())
  assert approx_equal(l.array_of_b(), g.array_of_b())
  assert approx_equal(l.c(), g.c())
  assert l.use_c()
  s = pickle.dumps(gaussian.sum((),()))
  l = pickle.loads(s)
  assert not l.use_c()
  exercise_gradient_dx(gaussian.sum(
    [5.5480], [10.4241], 0))
  exercise_gradient_dx(gaussian.sum(
   [2.657506,1.078079,1.490909,-4.241070,0.713791],
   [14.780758,0.776775,42.086842,-0.000294,0.239535],
   4.297983))
  exercise_integral_dx(gaussian.sum([5.5480], [10.4241]))
  exercise_integral_dx(gaussian.sum([5.5480], [10.4241], 3))
  exercise_integral_dx(gaussian.sum([5.5480], [0], 0))
  exercise_integral_dx(gaussian.sum([5.5480], [-0.01]))
  exercise_integral_dx(gaussian.sum(
    [2.657506,1.078079,1.490909,-4.241070,0.713791],
    [14.780758,0.776775,42.086842,-0.000294,0.239535],
    4.297983))
  g = gaussian.sum((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  s = StringIO()
  g.show(s)
  assert len(s.getvalue().split()) == 14
  g = gaussian.sum((3,-2,1,-4,5), (-.3,.2,-.1,.4,-.5))
  s = StringIO()
  g.show(s)
  assert len(s.getvalue().split()) == 12
  assert isinstance(g.sort(), gaussian.sum)
  assert approx_equal(g.sort().array_of_a(), (5,-4,3,-2,1))
  assert approx_equal(g.sort().array_of_b(), (-.5,.4,-.3,.2,-.1))
  assert not g.sort().use_c()
  g = gaussian.sum((1,2),(3,4),5)
  assert approx_equal(g.sort().array_of_a(), (2,1))
  assert approx_equal(g.sort().array_of_b(), (4,3))
  assert approx_equal(g.sort().c(), 5)
  assert g.sort().use_c()

def fit_finite_diff_gradients(gfit, x, eps=1.e-2):
  gr = flex.double()
  c = gfit.c()
  use_c = gfit.use_c()
  for i in range(gfit.n_terms()):
    t = []
    for seps in (eps, -eps):
      a = list(gfit.array_of_a())
      a[i] += seps
      t.append(
        gaussian.sum(a, gfit.array_of_b(), c, use_c).at_x(x))
    gr.append((t[0]-t[1])/(2*eps))
    t = []
    for seps in (eps, -eps):
      b = list(gfit.array_of_b())
      b[i] += seps
      t.append(
        gaussian.sum(gfit.array_of_a(), b, c, use_c).at_x(x))
    gr.append((t[0]-t[1])/(2*eps))
  if (use_c):
    t = []
    for seps in (eps, -eps):
      t.append(
        gaussian.sum(
          gfit.array_of_a(), gfit.array_of_b(), c+seps, use_c).at_x(x))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def fit_finite_diff_target_gradients(gfit, power, use_sigmas, eps=1.e-2):
  assert gfit.table_x().size() == 1
  weight = 1/gfit.table_sigmas()[0]**2
  gr = flex.double()
  c = gfit.c()
  use_c = gfit.use_c()
  for i in range(gfit.n_terms()):
    t = []
    for seps in (eps, -eps):
      a = list(gfit.array_of_a())
      a[i] += seps
      gf = gaussian.fit(
        gfit.table_x(),
        gfit.table_y(),
        gfit.table_sigmas(),
        gaussian.sum(a, gfit.array_of_b(), c, use_c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
    t = []
    for seps in (eps, -eps):
      b = list(gfit.array_of_b())
      b[i] += seps
      gf = gaussian.fit(
        gfit.table_x(),
        gfit.table_y(),
        gfit.table_sigmas(),
        gaussian.sum(gfit.array_of_a(), b, c, use_c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  if (use_c):
    t = []
    for seps in (eps, -eps):
      gf = gaussian.fit(
        gfit.table_x(),
        gfit.table_y(),
        gfit.table_sigmas(),
        gaussian.sum(gfit.array_of_a(), gfit.array_of_b(), c+seps, use_c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def exercise_fit():
  x = flex.double((0.1, 0.2, 0.5))
  y = flex.double((3,2,1))
  sigmas = flex.double((0.04,0.02,0.01))
  gf = gaussian.fit(
    x, y, sigmas,
    gaussian.sum((1,2), (4,5)))
  assert approx_equal(gf.array_of_a(), (1,2))
  assert approx_equal(gf.array_of_b(), (4,5))
  assert approx_equal(gf.c(), 0)
  assert not gf.use_c()
  assert approx_equal(gf.table_x(), x)
  assert approx_equal(gf.table_y(), y)
  assert approx_equal(gf.table_sigmas(), sigmas)
  assert approx_equal(gf.fitted_values(),
    [2.8632482881537511, 2.4896052951221748, 0.94088903489182252])
  reference_gaussian = gaussian.sum((1,2,3), (4,5,6))
  gf = gaussian.fit(
    x, reference_gaussian, sigmas,
    gaussian.sum((1,2), (4,5)))
  assert approx_equal(gf.array_of_a(), (1,2))
  assert approx_equal(gf.array_of_b(), (4,5))
  assert approx_equal(gf.c(), 0)
  assert approx_equal(gf.table_x(), x)
  assert approx_equal(gf.table_y(), reference_gaussian.at_x(x))
  assert approx_equal(gf.table_sigmas(), sigmas)
  assert isinstance(gf.sort(), gaussian.fit)
  assert gf.sort().table_x() == gf.table_x()
  assert gf.sort().table_y() == gf.table_y()
  assert gf.sort().table_sigmas() == gf.table_sigmas()
  assert approx_equal(gf.differences(), gf.at_x(x)-reference_gaussian.at_x(x))
  c_fit = gaussian.fit(
    flex.double([0.0, 0.066666666666666666, 0.13333333333333333,
                 0.2, 0.26666666666666666]),
    gaussian.sum(
      (2.657506, 1.078079, 1.490909, -4.2410698, 0.71379101),
      (14.780758, 0.776775, 42.086842, -0.000294, 0.239535),
      4.2979832),
    flex.double(5, 0.0005),
    gaussian.sum(
      (1.1423916, 4.1728425, 0.61716694),
      (0.50733125, 14.002512, 41.978928)))
  differences = flex.double([-0.064797341823577881, 0.003608505180995536,
    0.098159179757290715, 0.060724224581695019, -0.10766283796372011])
  assert approx_equal(c_fit.differences(), differences)
  assert approx_equal(c_fit.significant_relative_errors(),
    [0.0107212, 0.0005581, 0.0213236, 0.0169304, 0.0385142])
  gf = gaussian.fit(
    x, reference_gaussian, flex.double(x.size(), 1),
    gaussian.sum((1,2), (4,5)))
  assert list(gf.bound_flags(False, False)) == [False,False,False,False]
  assert list(gf.bound_flags(True, False)) == [True,False,True,False]
  assert list(gf.bound_flags(False, True)) == [False,True,False,True]
  sgf = gf.apply_shifts(flex.double((3,-3,4,6)), True)
  assert approx_equal(sgf.array_of_a(), (1+3,2+4))
  assert approx_equal(sgf.array_of_b(),
    ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgf.c(), 0)
  assert not sgf.use_c()
  sgf = gf.apply_shifts(flex.double((3,-3,4,6)), False)
  assert approx_equal(sgf.array_of_a(), (1+3,2+4))
  assert approx_equal(sgf.array_of_b(), (4-3,5+6))
  assert approx_equal(sgf.c(), 0)
  assert not sgf.use_c()
  differences = sgf.differences()
  for use_sigmas in [False, True]:
    assert approx_equal(sgf.target_function(2, use_sigmas, differences),
      25.0320634)
    assert approx_equal(sgf.target_function(4, use_sigmas, differences),
      256.2682575)
    assert approx_equal(
      sgf.gradients_d_abc(2, use_sigmas, differences),
      [15.6539271, -4.1090114, 10.4562306, -1.6376781])
  gfc = gaussian.fit(
    x, reference_gaussian, flex.double(x.size(), 1),
    gaussian.sum((1,2), (4,5), 6))
  assert list(gfc.bound_flags(False, False)) == [False,False,False,False,False]
  assert list(gfc.bound_flags(True, False)) == [True,False,True,False,True]
  assert list(gfc.bound_flags(False, True)) == [False,True,False,True,False]
  sgfc = gfc.apply_shifts(flex.double((3,-3,4,6,-5)), True)
  assert approx_equal(sgfc.array_of_a(), (1+3,2+4))
  assert approx_equal(sgfc.array_of_b(),
    ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgfc.c(), 6-5)
  assert sgfc.use_c()
  sgfc = gfc.apply_shifts(flex.double((3,-3,4,6,-5)), False)
  assert approx_equal(sgfc.array_of_a(), (1+3,2+4))
  assert approx_equal(sgfc.array_of_b(), (4-3,5+6))
  assert approx_equal(sgfc.c(), 6-5)
  assert sgfc.use_c()
  differences = sgfc.differences()
  for use_sigmas in [False, True]:
    assert approx_equal(sgfc.target_function(2, use_sigmas, differences),
      44.8181444)
    assert approx_equal(sgfc.target_function(4, use_sigmas, differences),
      757.3160329)
    assert approx_equal(
      sgfc.gradients_d_abc(2, use_sigmas, differences),
      [21.1132071, -6.0532695, 13.6638274, -2.2460994, 22.7860809])
  differences = c_fit.differences()
  gabc = c_fit.gradients_d_abc(2, False, differences)
  assert approx_equal(
    gabc,
    [-0.016525391425206391, 0.0074465239375589107, 0.020055876723667564,
     0.00054794635257838251, -0.018754011379726425, -0.0011194004809549143])
  assert approx_equal(
    c_fit.gradients_d_shifts(flex.double((0.1,0.4,0.2,0.5,0.3,0.6)), gabc),
    [-0.0165254, 0.01656512, 0.0200559, 0.0046488, -0.0187540, -0.0158487])
  g5c = gaussian.sum(
    (2.657505989074707, 1.0780789852142334, 1.4909089803695679,
     -4.2410697937011719, 0.71379101276397705),
    (14.780757904052734, 0.77677500247955322, 42.086841583251953,
     -0.00029399999766610563, 0.23953500390052795),
    4.2979831695556641)
  for include_constant_term in (False, True):
    a = flex.double(g5c.array_of_a())
    b = flex.double(g5c.array_of_b())
    permutation = flex.sort_permutation(data=flex.abs(a), reverse=True)[:4]
    gf = gaussian.fit(
      flex.double([0]),
      g5c,
      flex.double(1, 1),
      gaussian.sum(
        iter(a.select(permutation)),
        iter(b.select(permutation)), 0, include_constant_term))
    assert approx_equal(gf.differences(), [-5.01177418232])
    shifts = flex.double(8,-1)
    if (include_constant_term): shifts.append(-.2)
    sgf = gf.apply_shifts(shifts, False)
    assert approx_equal(sgf.array_of_a(),
                        [-5.2410698, 1.657506, 0.49090898, 0.078078985])
    assert approx_equal(sgf.array_of_b(),
                        [-1.0002940, 13.780758, 41.086842, -0.223225])
    if (include_constant_term):
      assert approx_equal(sgf.c(), -.2)
    expected_gradients = [1,0,1,0,1,0,1,0]
    if (include_constant_term): expected_gradients.append(1)
    assert approx_equal(
      fit_finite_diff_gradients(sgf, 0),
      expected_gradients,
      eps=1.e-4)
    for i in range(10):
      gf = gaussian.fit(
        flex.double([i / 10.]),
        g5c,
        flex.double(1, 1),
        sgf)
      differences = flex.double([0.5])
      assert approx_equal(
        gf.gradients_d_abc(2, False, differences),
        fit_finite_diff_gradients(gf, gf.table_x()[0]),
        eps=1.e-3)
      for sigma in [0.04,0.02,0.01]:
        gf = gaussian.fit(
          flex.double([i / 20.]),
          g5c,
          flex.double([sigma]),
          sgf)
        for power in [2,4]:
          for use_sigmas in [False, True]:
            differences = gf.differences()
            an=gf.gradients_d_abc(power, use_sigmas, differences)
            fi=fit_finite_diff_target_gradients(gf, power, use_sigmas)
            assert eps_eq(an, fi, eps=1.e-3)

carbon_s_y_table = [
0.00, 6.000, 0.01, 5.990, 0.02, 5.958, 0.03, 5.907, 0.04, 5.837, 0.05, 5.749,
0.06, 5.645, 0.07, 5.526, 0.08, 5.396, 0.09, 5.255, 0.10, 5.107, 0.11, 4.952,
0.12, 4.794, 0.13, 4.633, 0.14, 4.472, 0.15, 4.311, 0.16, 4.153, 0.17, 3.998,
0.18, 3.847, 0.19, 3.701, 0.20, 3.560, 0.22, 3.297, 0.24, 3.058, 0.25, 2.949,
0.26, 2.846, 0.28, 2.658, 0.30, 2.494, 0.32, 2.351, 0.34, 2.227, 0.35, 2.171,
0.36, 2.120, 0.38, 2.028, 0.40, 1.948, 0.42, 1.880, 0.44, 1.821, 0.45, 1.794,
0.46, 1.770, 0.48, 1.725, 0.50, 1.685, 0.55, 1.603, 0.60, 1.537, 0.65, 1.479,
0.70, 1.426, 0.80, 1.322, 0.90, 1.219, 1.00, 1.114, 1.10, 1.012, 1.20, 0.914,
1.30, 0.822, 1.40, 0.736, 1.50, 0.659, 1.60, 0.588, 1.70, 0.525, 1.80, 0.468,
1.90, 0.418, 2.00, 0.373, 2.50, 0.216, 3.00, 0.130, 3.50, 0.081, 4.00, 0.053,
5.00, 0.025, 6.00, 0.013]

class tabulated_fit:

  def __init__(self, limit, coefficients):
    self.limit = limit
    self.coefficients = coefficients

carbon_fit_6 = tabulated_fit(6.0, [
2.18188567686, 13.4533708328,
1.77612377639, 32.5790123523,
1.08772011297, 0.747293264573,
0.641460989931, 0.251251498175,
0.207885994451, 80.9799313275,
0.105219184507, 0.0587297979816])

carbon_fit_5 = tabulated_fit(6.0, [
2.65463431663, 14.7665037505,
1.49420264709, 42.0409767208,
1.05563210943, 0.780856499884,
0.688021531597, 0.258963998784,
0.104681246572, 0.0579465611728])

carbon_fit_4 = tabulated_fit(3.0, [
2.21557580709, 12.7523000206,
1.98306066831, 36.4905110196,
1.31636728472, 0.632825354093,
0.480812064621, 0.148079120135])

carbon_fit_3 = tabulated_fit(1.4, [
2.51340127252, 31.8053433708,
1.74867019409, 0.445605499982,
1.72398202356, 10.5831679451])

carbon_fit_2 = tabulated_fit(0.5, [
3.54355550695, 25.6239838191,
2.42579673735, 1.50364460774])

carbon_fit_1 = tabulated_fit(0.15, [
5.96792806111, 14.8957682987])

carbon_it1992 = tabulated_fit(2.0, [
2.31000, 20.8439,
1.02000, 10.2075,
1.58860, 0.568700,
0.865000, 51.6512,
0.215600])

carbon_wk1995 = tabulated_fit(6.0, [
2.657506, 14.780758,
1.078079, 0.776775,
1.490909, 42.086842,
-4.241070, -0.000294,
0.713791, 0.239535,
4.297983])

class carbon_fit(immoptibox_ports.test_function):

  def __init__(self, tab_fit, perturb, verbose):
    self.tab_fit = tab_fit
    self.perturb = perturb
    self.verbose = verbose
    carbon_ss = flex.double(carbon_s_y_table)[0::2]
    carbon_ys = flex.double(carbon_s_y_table)[1::2]
    selection = carbon_ss <= tab_fit.limit + 1.e-3
    self.fit = gaussian.fit(
      carbon_ss.select(selection),
      carbon_ys.select(selection),
      flex.double(selection.count(True), 1),
      gaussian.sum(flex.double(tab_fit.coefficients)))
    n = self.fit.n_parameters()
    immoptibox_ports.test_function.__init__(self,
      m=self.fit.table_x().size(),
      n=n,
      check_with_finite_differences=(n <= 6 or n == 9),
      verbose=verbose)

  def initialization(self):
    self.x0 = self.fit.parameters()
    self.capital_f_x_star = 0.5*self.f(x=self.x0).norm()**2
    if (self.perturb):
      mersenne_twister = flex.mersenne_twister(seed=0)
      self.x0 *= 1 + mersenne_twister.random_double(
        size=self.x0.size(), factor=0.01)
    self.tau0 = 1e-8
    self.delta0 = 10
    self.x_star = None

  def label(self):
    return "carbon_fit(n=%d, perturb=%s)" % (
      self.fit.n_parameters(), str(self.perturb))

  def check_minimized_capital_f_x_star(self, f_x_star, tolerance=1.e-3):
    capital_f_x_star = 0.5*f_x_star.norm()**2
    if (capital_f_x_star > self.capital_f_x_star):
      assert capital_f_x_star < tolerance, (
        capital_f_x_star, self.capital_f_x_star)
      if (self.verbose):
        print("  WARNING: minimization converged to larger residual", \
          "than original solution:")
        print("    original:", self.capital_f_x_star)
      assert self.perturb

  def f(self, x):
    fit = gaussian.fit(
      self.fit.table_x(), self.fit.table_y(), self.fit.table_sigmas(),
      gaussian.sum(x))
    return fit.differences()

  def jacobian_analytical(self, x):
    fit = gaussian.fit(
      self.fit.table_x(), self.fit.table_y(), self.fit.table_sigmas(),
      gaussian.sum(x))
    return fit.least_squares_jacobian_abc()

  def hessian_analytical(self, x):
    j = self.jacobian_analytical(x=x)
    fit = gaussian.fit(
      self.fit.table_x(), self.fit.table_y(), self.fit.table_sigmas(),
      gaussian.sum(x))
    return fit.least_squares_hessian_abc_as_packed_u() \
      .matrix_packed_u_as_symmetric()

def exercise_fit_jacobian_and_hessian(verbose):
  for tab_fit in [carbon_fit_1, carbon_fit_2, carbon_fit_3,
                  carbon_fit_4, carbon_fit_5, carbon_fit_6,
                  carbon_it1992, carbon_wk1995]:
    for perturb in [False, True]:
      carbon_fit(tab_fit=tab_fit, perturb=perturb, verbose=verbose)

def run():
  exercise_term()
  exercise_sum()
  exercise_fit()
  exercise_fit_jacobian_and_hessian(verbose="--verbose" in sys.argv[1:])
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
