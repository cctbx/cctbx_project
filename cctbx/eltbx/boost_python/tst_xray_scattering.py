from cctbx.eltbx import xray_scattering
from scitbx.array_family import flex
from scitbx.test_utils import approx_equal
import StringIO
import pickle
import string

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
    [-1, 13.780758, 41.086842, -0.223225]).gradients_at_d_star_sq(0.81)
  assert approx_equal(grg.a(),
    [1.2244601, 0.06138416, 0.00024357505, 1.0462403])
  assert approx_equal(grg.b(),
    [1.2995399, -0.020603284, -2.4213568e-05, -0.016542099])
  assert approx_equal(grg.c(), 1)
  s = StringIO.StringIO()
  grg.show(s)
  assert len(s.getvalue().split()) == 12

def finite_diff_gradients(diff_gaussian, d_star_sq, eps=1.e-2):
  gr = flex.double()
  for i in xrange(diff_gaussian.n_ab()):
    t = []
    for seps in (eps, -eps):
      a = list(diff_gaussian.a())
      a[i] += seps
      t.append(
        xray_scattering.gaussian(a, diff_gaussian.b()).at_d_star_sq(d_star_sq))
    gr.append((t[0]-t[1])/(2*eps))
  for i in xrange(diff_gaussian.n_ab()):
    t = []
    for seps in (eps, -eps):
      b = list(diff_gaussian.b())
      b[i] += seps
      t.append(
        xray_scattering.gaussian(diff_gaussian.a(), b).at_d_star_sq(d_star_sq))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def finite_diff_target_gradients(diff_gaussian, d_star_sq, weight, eps=1.e-2):
  gr = flex.double()
  for i in xrange(diff_gaussian.n_ab()):
    t = []
    for seps in (eps, -eps):
      a = list(diff_gaussian.a())
      a[i] += seps
      t.append(xray_scattering.difference_gaussian(
        diff_gaussian.reference_gaussian(),
        xray_scattering.gaussian(
          a, diff_gaussian.b())).target_at_d_star_sq(d_star_sq, weight))
    gr.append((t[0]-t[1])/(2*eps))
  for i in xrange(diff_gaussian.n_ab()):
    t = []
    for seps in (eps, -eps):
      b = list(diff_gaussian.b())
      b[i] += seps
      t.append(xray_scattering.difference_gaussian(
        diff_gaussian.reference_gaussian(),
        xray_scattering.gaussian(
          diff_gaussian.a(), b)).target_at_d_star_sq(d_star_sq, weight))
    gr.append((t[0]-t[1])/(2*eps))
  return gr

def exercise_difference_gaussian():
  dg = xray_scattering.difference_gaussian(
    xray_scattering.gaussian((1,2,3), (4,5,6), 7),
    xray_scattering.gaussian((1,2), (4,5)))
  assert approx_equal(dg.a(), (1,2))
  assert approx_equal(dg.b(), (4,5))
  assert approx_equal(dg.c(), 0)
  assert approx_equal(dg.reference_gaussian().a(), (1,2,3))
  assert approx_equal(dg.reference_gaussian().b(), (4,5,6))
  assert approx_equal(dg.reference_gaussian().c(), 7)
  sdg = dg.apply_shifts(flex.double((3,4,-5,6)), 1)
  assert approx_equal(sdg.a(), (1+3,2+4))
  assert approx_equal(sdg.b(), (1,5+6))
  assert approx_equal(sdg.c(), 0)
  assert approx_equal(sdg.target_term_at_d_star_sq(0), -3)
  assert approx_equal(sdg.target_at_d_star_sq(0, 1), 9)
  d_star_sq = flex.double((0,.5,1))
  weights = flex.double((1,1,1))
  target_terms = sdg.target_terms_at_points(d_star_sq)
  assert approx_equal(target_terms, [-3.0, -5.0471280, -5.1115092])
  assert approx_equal(
    sdg.sum_of_gradients_at_points(d_star_sq, weights, target_terms),
    [-22.8698443, -9.2057633, 12.4157696, 2.8944743])
  dg = xray_scattering.difference_gaussian(
    xray_scattering.gaussian(
      (2.657506, 1.078079, 1.490909, -4.2410698, 0.71379101),
      (14.780758, 0.776775, 42.086842, -0.000294, 0.239535),
      4.2979832),
    xray_scattering.gaussian(
      (1.1423916, 4.1728425, 0.61716694),
      (0.50733125, 14.002512, 41.978928)))
  d_star_sq = flex.double([0.0, 0.017777777777777778, 0.071111111111111111,
    0.16000000000000003, 0.28444444444444444])
  weights = flex.double(5, 1)
  target_terms = flex.double([-0.064797341823577881, 0.003608505180995536,
    0.098159179757290715, 0.060724224581695019, -0.10766283796372011])
  assert approx_equal(dg.target_terms_at_points(d_star_sq), target_terms)
  assert approx_equal(
    dg.sum_of_gradients_at_points(d_star_sq, weights, target_terms),
    [-0.016525391425206391, 0.020055876723667564, -0.018754011379726425,
    0.0074465239375589107, 0.00054794635257838251, -0.0011194004809549143])
  g5c = xray_scattering.wk1995("C")
  a = flex.double(g5c.a())
  b = flex.double(g5c.b())
  permutation = flex.sort_permutation(flex.abs(a), 1)[:4]
  gdiff = xray_scattering.difference_gaussian(
    g5c.fetch(),
    xray_scattering.gaussian(
      iter(a.select(permutation)),
      iter(b.select(permutation))))
  assert approx_equal(gdiff.target_term_at_d_star_sq(0), -5.01177418232)
  assert approx_equal(gdiff.target_at_d_star_sq(0, 1), 25.1178804546)
  gshifted = gdiff.apply_shifts(flex.double(8,-1), -1)
  assert approx_equal(gshifted.a(),
                      [-5.2410698, 1.657506, 0.49090898, 0.078078985])
  assert approx_equal(gshifted.b(),
                      [-1, 13.780758, 41.086842, -0.223225])
  assert approx_equal(finite_diff_gradients(gshifted, 0),
                      [1,1,1,1,0,0,0,0], eps=1.e-4)
  for weight in [1,2,3]:
    for i in xrange(10):
      d_star_sq = flex.double([(i / 10.)**2])
      weights = flex.double([weight])
      tt = flex.double([0.5])
      assert approx_equal(
        gshifted.sum_of_gradients_at_points(d_star_sq, flex.double(1,1), tt),
        finite_diff_gradients(gshifted, d_star_sq[0]), eps=1.e-4)
      tt = gshifted.target_terms_at_points(d_star_sq)
      assert approx_equal(
        gshifted.sum_of_gradients_at_points(d_star_sq, weights, tt),
        finite_diff_target_gradients(gshifted, d_star_sq[0], weights[0]),
        eps=1.e-2)

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
  exercise_difference_gaussian()
  exercise_it1992()
  exercise_wk1995()
  ensure_common_symbols()
  ensure_correct_element_symbol()
  print "OK"

if (__name__ == "__main__"):
  run()
