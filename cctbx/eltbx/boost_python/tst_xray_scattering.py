from cctbx.eltbx import xray_scattering
from scitbx.test_utils import approx_equal
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
    [5.5480], [10.4241], 0))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [10.4241], 3))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [0], 0))
  exercise_gaussian_integral(xray_scattering.gaussian(
    [5.5480], [-0.01], 0))
  exercise_gaussian_integral(xray_scattering.gaussian(
   [2.657506,1.078079,1.490909,-4.241070,0.713791],
   [14.780758,0.776775,42.086842,-0.000294,0.239535],
   4.297983))

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
  exercise_it1992()
  exercise_wk1995()
  ensure_common_symbols()
  ensure_correct_element_symbol()
  print "OK"

if (__name__ == "__main__"):
  run()
