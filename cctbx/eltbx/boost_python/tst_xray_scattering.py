from cctbx.eltbx import xray_scattering
from scitbx.array_family import flex
from scitbx.test_utils import approx_equal, eps_eq
import StringIO
import pickle
import string
import math

def gaussian_finite_gradient(gaussian, d_star, eps=1.e-6):
  if (d_star == 0): return 0
  assert d_star >= eps
  tm = gaussian.at_d_star_sq((d_star-eps)**2)
  tp = gaussian.at_d_star_sq((d_star+eps)**2)
  return (tp-tm)/(2*eps)

def exercise_gaussian_gradient(gaussian, d_min=0.5, n_points=50):
  d_star_max = 1./d_min
  for i in xrange(n_points+1):
    d_star = d_star_max * i / n_points
    grad_finite = gaussian_finite_gradient(gaussian, d_star)
    grad_analytical = gaussian.gradient_at_d_star(d_star)
    if (0):
      print d_star, grad_finite, grad_analytical
    assert abs(grad_finite - grad_analytical) < 1.e-6

def exercise_gaussian_integral(gaussian, d_min=0.5, n_points=1000):
  d_star_max = 1./d_min
  d_step = d_star_max / n_points
  numerical_integral = 0
  for i in xrange(n_points+1):
    d_star = d_star_max * i / n_points
    new_value = gaussian.at_d_star_sq(d_star**2)
    if (i):
      numerical_integral += (prev_value + new_value) * .5
    prev_value = new_value
    analytical_integral = gaussian.integral_at_d_star(d_star)
    if (0):
      print d_star, numerical_integral*d_step, analytical_integral
    assert abs(numerical_integral*d_step - analytical_integral) < 1.e-4

def exercise_gaussian():
  c = xray_scattering.gaussian(0)
  assert c.n_ab() == 0
  assert c.a() == ()
  assert c.b() == ()
  assert approx_equal(c.c(), 0)
  assert c.all_zero()
  c = xray_scattering.gaussian(1)
  assert c.n_ab() == 0
  assert c.a() == ()
  assert c.b() == ()
  assert approx_equal(c.c(), 1)
  assert not c.all_zero()
  c = xray_scattering.gaussian((), ())
  assert c.n_ab() == 0
  assert c.a() == ()
  assert c.b() == ()
  assert c.c() == 0
  c = xray_scattering.gaussian((), (), -2)
  assert c.n_ab() == 0
  assert c.a() == ()
  assert c.b() == ()
  assert approx_equal(c.c(), -2)
  c = xray_scattering.gaussian((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  assert c.n_ab() == 5
  assert approx_equal(c.a(),(1,-2,3,-4,5))
  assert approx_equal(c.b(),(-.1,.2,-.3,.4,-.5))
  assert approx_equal(c.c(), 6)
  assert approx_equal(c.at_d_star_sq(3), 13.4251206)
  assert approx_equal(c.at_d_star_sq(flex.double([2,3])),
    [11.8723031, 13.4251206])
  s = pickle.dumps(c)
  l = pickle.loads(s)
  assert l.n_ab() == c.n_ab()
  assert approx_equal(l.a(), c.a())
  assert approx_equal(l.b(), c.b())
  assert approx_equal(l.c(), c.c())
  exercise_gaussian_gradient(xray_scattering.gaussian(
    [5.5480], [10.4241], 0))
  exercise_gaussian_gradient(xray_scattering.gaussian(
   [2.657506,1.078079,1.490909,-4.241070,0.713791],
   [14.780758,0.776775,42.086842,-0.000294,0.239535],
   4.297983))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [10.4241]))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [10.4241], 3))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [0], 0))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [-0.01]))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [2.657506,1.078079,1.490909,-4.241070,0.713791],
    [14.780758,0.776775,42.086842,-0.000294,0.239535],
    4.297983))
  grg = xray_scattering.gaussian(
    [-5.2410698, 1.657506, 0.49090898, 0.078078985],
    [-1, 13.780758, 41.086842, -0.223225]).gradients_at_stol(0.45)
  assert approx_equal(grg.a(),
    [1.2244601, 0.06138416, 0.00024357505, 1.0462403])
  assert approx_equal(grg.b(),
    [1.2995399, -0.020603284, -2.4213568e-05, -0.016542099])
  assert approx_equal(grg.c(), 1)
  s = StringIO.StringIO()
  grg.show(s)
  assert len(s.getvalue().split()) == 12

def finite_diff_gradients(diff_gaussian, include_constant_term, stol,
                          eps=1.e-2):
  gr = flex.double()
  c = diff_gaussian.c()
  for i in xrange(diff_gaussian.n_ab()):
    t = []
    for seps in (eps, -eps):
      a = list(diff_gaussian.a())
      a[i] += seps
      t.append(
        xray_scattering.gaussian(a, diff_gaussian.b(), c).at_stol(stol))
    gr.append((t[0]-t[1])/(2*eps))
  for i in xrange(diff_gaussian.n_ab()):
    t = []
    for seps in (eps, -eps):
      b = list(diff_gaussian.b())
      b[i] += seps
      t.append(
        xray_scattering.gaussian(diff_gaussian.a(), b, c).at_stol(stol))
    gr.append((t[0]-t[1])/(2*eps))
  if (include_constant_term):
    t = []
    for seps in (eps, -eps):
      t.append(
        xray_scattering.gaussian(
          diff_gaussian.a(), diff_gaussian.b(), c+seps).at_stol(stol))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def finite_diff_target_gradients(gaussian_fit, include_constant_term,
                                 power, use_sigmas, eps=1.e-2):
  assert gaussian_fit.stols().size() == 1
  weight = 1/gaussian_fit.sigmas()[0]**2
  gr = flex.double()
  c = gaussian_fit.c()
  for i in xrange(gaussian_fit.n_ab()):
    t = []
    for seps in (eps, -eps):
      a = list(gaussian_fit.a())
      a[i] += seps
      gf = xray_scattering.gaussian_fit(
        gaussian_fit.stols(),
        gaussian_fit.target_values(),
        gaussian_fit.sigmas(),
        xray_scattering.gaussian(a, gaussian_fit.b(), c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  for i in xrange(gaussian_fit.n_ab()):
    t = []
    for seps in (eps, -eps):
      b = list(gaussian_fit.b())
      b[i] += seps
      gf = xray_scattering.gaussian_fit(
        gaussian_fit.stols(),
        gaussian_fit.target_values(),
        gaussian_fit.sigmas(),
        xray_scattering.gaussian(
          gaussian_fit.a(), b, c))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  if (include_constant_term):
    t = []
    for seps in (eps, -eps):
      gf = xray_scattering.gaussian_fit(
        gaussian_fit.stols(),
        gaussian_fit.target_values(),
        gaussian_fit.sigmas(),
        xray_scattering.gaussian(
          gaussian_fit.a(), gaussian_fit.b(), c+seps))
      t.append(gf.target_function(power, use_sigmas, gf.differences()))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def exercise_gaussian_fit():
  stols = flex.double((0.1, 0.2, 0.5))
  values = flex.double((3,2,1))
  sigmas = flex.double((0.04,0.02,0.01))
  gf = xray_scattering.gaussian_fit(
    stols, values, sigmas,
    xray_scattering.gaussian((1,2), (4,5)))
  assert approx_equal(gf.a(), (1,2))
  assert approx_equal(gf.b(), (4,5))
  assert approx_equal(gf.c(), 0)
  assert approx_equal(gf.stols(), stols)
  assert approx_equal(gf.target_values(), values)
  assert approx_equal(gf.sigmas(), sigmas)
  assert approx_equal(gf.fitted_values(),
    [2.8632482881537511, 2.4896052951221748, 0.94088903489182252])
  gf = xray_scattering.gaussian_fit(
    stols, values, flex.double(),
    xray_scattering.gaussian((1,2), (4,5)))
  assert approx_equal(gf.sigmas(), [1,1,1])
  reference_gaussian = xray_scattering.gaussian((1,2,3), (4,5,6))
  gf = xray_scattering.gaussian_fit(
    stols, reference_gaussian, sigmas,
    xray_scattering.gaussian((1,2), (4,5)))
  assert approx_equal(gf.a(), (1,2))
  assert approx_equal(gf.b(), (4,5))
  assert approx_equal(gf.c(), 0)
  assert approx_equal(gf.stols(), stols)
  assert approx_equal(gf.target_values(),
    [reference_gaussian.at_stol(stol) for stol in stols])
  assert approx_equal(gf.sigmas(), sigmas)
  gf = xray_scattering.gaussian_fit(
    stols, reference_gaussian, flex.double(),
    xray_scattering.gaussian((1,2), (4,5)))
  assert approx_equal(gf.sigmas(), [1,1,1])
  sgf = gf.apply_shifts(flex.double((3,4,-3,6)), 0001)
  assert approx_equal(sgf.a(), (1+3,2+4))
  assert approx_equal(sgf.b(), ((math.sqrt(4)-3)**2,(math.sqrt(5)+6)**2))
  assert approx_equal(sgf.c(), 0)
  sgf = gf.apply_shifts(flex.double((3,4,-3,6)), 00000)
  assert approx_equal(sgf.a(), (1+3,2+4))
  assert approx_equal(sgf.b(), (4-3,5+6))
  assert approx_equal(sgf.c(), 0)
  differences = sgf.differences()
  assert approx_equal(differences,
    [sgf.at_stol(stol)-reference_gaussian.at_stol(stol) for stol in stols])
  for use_sigmas in [00000, 0001]:
    assert approx_equal(sgf.target_function(2, use_sigmas, differences),
      25.0320634)
    assert approx_equal(sgf.target_function(4, use_sigmas, differences),
      256.2682575)
    assert approx_equal(
      sgf.gradients_w_r_t_abc(2, use_sigmas, differences, 00000),
      [15.6539271, 10.4562306, -4.1090114, -1.6376781])
  gf = xray_scattering.gaussian_fit(
    flex.double([0.0, 0.066666666666666666, 0.13333333333333333,
                 0.2, 0.26666666666666666]),
    xray_scattering.gaussian(
      (2.657506, 1.078079, 1.490909, -4.2410698, 0.71379101),
      (14.780758, 0.776775, 42.086842, -0.000294, 0.239535),
      4.2979832),
    flex.double(),
    xray_scattering.gaussian(
      (1.1423916, 4.1728425, 0.61716694),
      (0.50733125, 14.002512, 41.978928)))
  differences = flex.double([-0.064797341823577881, 0.003608505180995536,
    0.098159179757290715, 0.060724224581695019, -0.10766283796372011])
  assert approx_equal(gf.differences(), differences)
  assert approx_equal(
    gf.gradients_w_r_t_abc(2, 00000, differences, 00000),
    [-0.016525391425206391, 0.020055876723667564, -0.018754011379726425,
    0.0074465239375589107, 0.00054794635257838251, -0.0011194004809549143])
  for include_constant_term in (00000, 0001):
    g5c = xray_scattering.wk1995("C")
    a = flex.double(g5c.a())
    b = flex.double(g5c.b())
    permutation = flex.sort_permutation(flex.abs(a), 1)[:4]
    gf = xray_scattering.gaussian_fit(
      flex.double([0]),
      g5c.fetch(),
      flex.double(),
      xray_scattering.gaussian(
        iter(a.select(permutation)),
        iter(b.select(permutation)), 0))
    assert approx_equal(gf.differences(), [-5.01177418232])
    shifts = flex.double(8,-1)
    if (include_constant_term): shifts.append(-.2)
    sgf = gf.apply_shifts(shifts, 00000)
    assert approx_equal(sgf.a(),
                        [-5.2410698, 1.657506, 0.49090898, 0.078078985])
    assert approx_equal(sgf.b(),
                        [-1.0002940, 13.780758, 41.086842, -0.223225])
    if (include_constant_term):
      assert approx_equal(sgf.c(), -.2)
    expected_gradients = [1,1,1,1,0,0,0,0]
    if (include_constant_term): expected_gradients.append(1)
    assert approx_equal(
      finite_diff_gradients(sgf, include_constant_term, 0),
      expected_gradients,
      eps=1.e-4)
    for i in xrange(10):
      gf = xray_scattering.gaussian_fit(
        flex.double([i / 10.]),
        g5c.fetch(),
        flex.double(),
        sgf)
      differences = flex.double([0.5])
      assert approx_equal(
        gf.gradients_w_r_t_abc(2, 00000, differences, include_constant_term),
        finite_diff_gradients(gf, include_constant_term, gf.stols()[0]),
        eps=1.e-4)
      for sigma in [0.04,0.02,0.01]:
        gf = xray_scattering.gaussian_fit(
          flex.double([i / 20.]),
          g5c.fetch(),
          flex.double([sigma]),
          sgf)
        for power in [2,4]:
          for use_sigmas in [00000, 0001]:
            differences = gf.differences()
            an = gf.gradients_w_r_t_abc(
              power, use_sigmas, differences, include_constant_term)
            fi = finite_diff_target_gradients(
              gf, include_constant_term, power, use_sigmas)
            assert eps_eq(an, fi, eps=1.e-3)

def exercise_it1992():
  c = xray_scattering.it1992("c1")
  assert c.label() == "C"
  assert approx_equal(c.a(), (2.31000, 1.02000, 1.58860, 0.865000))
  assert approx_equal(c.b(), (20.8439, 10.2075, 0.568700, 51.6512))
  assert approx_equal(c.c(), 0.215600)
  assert approx_equal(c.at_stol_sq(0), 5.99919997156)
  assert approx_equal(c.at_stol_sq(1./9), 2.26575563201)
  assert approx_equal(c.at_stol(1./9), 4.93537567523)
  assert approx_equal(c.at_d_star_sq(1./9), 4.04815863088)
  c = xray_scattering.it1992("yb2+", 1)
  assert c.label() == "Yb2+"
  assert approx_equal(c.a()[0], 28.1209)
  assert approx_equal(c.b()[3], 20.3900)
  assert approx_equal(c.c(), 3.70983)
  c = xray_scattering.it1992("  YB3+")
  assert c.label() == "Yb3+"
  assert approx_equal(c.a()[0], 27.8917)
  n = 0
  for c in xray_scattering.it1992_iterator():
    n += 1
    if (n == 215):
      assert c.label() == "Cf"
    d = xray_scattering.it1992(c.label(), 1)
    assert d.label() == c.label()
  assert n == 215
  i = xray_scattering.it1992_iterator()
  j = iter(i)
  assert i is j
  t = xray_scattering.it1992("Si")
  c = t.fetch()
  assert c.n_ab() == 4
  assert c.a() == t.a()
  assert c.b() == t.b()
  assert c.c() == t.c()

def exercise_wk1995():
  c = xray_scattering.wk1995("c1")
  assert c.label() == "C"
  assert approx_equal(c.a(), (2.657506,1.078079,1.490909,-4.241070,0.713791))
  assert approx_equal(c.b(), (14.780758,0.776775,42.086842,-0.000294,0.239535))
  assert approx_equal(c.c(), 4.297983)
  assert approx_equal(c.at_stol_sq(0), 5.99719834328)
  assert approx_equal(c.at_stol_sq(1./9), 2.26895371584)
  assert approx_equal(c.at_stol(1./9), 4.93735084739)
  assert approx_equal(c.at_d_star_sq(1./9), 4.04679561237)
  c = xray_scattering.wk1995("yb2+", 1)
  assert c.label() == "Yb2+"
  assert approx_equal(c.a()[0], 28.443794)
  assert approx_equal(c.b()[4], 0.001463)
  assert approx_equal(c.c(), -23.214935)
  c = xray_scattering.wk1995("  YB3+")
  assert c.label() == "Yb3+"
  assert approx_equal(c.a()[0], 28.191629)
  n = 0
  for c in xray_scattering.wk1995_iterator():
    n += 1
    if (n == 215):
      assert c.label() == "Pu6+"
    d = xray_scattering.wk1995(c.label(), 1)
    assert d.label() == c.label()
  assert n == 215, n
  i = xray_scattering.wk1995_iterator()
  j = iter(i)
  assert i is j
  t = xray_scattering.wk1995("Si")
  c = t.fetch()
  assert c.n_ab() == 5
  assert c.a() == t.a()
  assert c.b() == t.b()
  assert c.c() == t.c()

def ensure_common_symbols():
  lbl_it = []
  for c in xray_scattering.it1992_iterator(): lbl_it.append(c.label())
  lbl_it.sort()
  lbl_wk = []
  for c in xray_scattering.wk1995_iterator(): lbl_wk.append(c.label())
  lbl_wk.sort()
  assert lbl_it == lbl_wk

def ensure_correct_element_symbol():
  from cctbx.eltbx import tiny_pse
  for g in xray_scattering.it1992_iterator():
    s = g.label()
    if (s == "Cval"):
      e = "C"
    else:
      e = s[:2]
      if (len(e) > 1 and e[1] not in string.letters):
        e = e[:1]
    assert tiny_pse.table(s).symbol() == e
    assert tiny_pse.table(s.lower()).symbol() == e
    assert tiny_pse.table(s.upper()).symbol() == e

def run():
  exercise_gaussian()
  exercise_gaussian_fit()
  exercise_it1992()
  exercise_wk1995()
  ensure_common_symbols()
  ensure_correct_element_symbol()
  print "OK"

if (__name__ == "__main__"):
  run()
