from cctbx.eltbx import xray_scattering
from scitbx.test_utils import approx_equal
import pickle
import math
import string

def exercise_gaussian():
  g = xray_scattering.gaussian(0)
  assert g.n_terms() == 0
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  assert g.n_parameters() == 1
  g = xray_scattering.gaussian(0, 00000)
  assert g.n_terms() == 0
  assert approx_equal(g.c(), 0)
  assert not g.use_c()
  assert g.n_parameters() == 0
  g = xray_scattering.gaussian(1)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), 1)
  assert g.n_parameters() == 1
  g = xray_scattering.gaussian((), ())
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert g.c() == 0
  g = xray_scattering.gaussian((), (), -2)
  assert g.n_terms() == 0
  assert g.array_of_a() == ()
  assert g.array_of_b() == ()
  assert approx_equal(g.c(), -2)
  g = xray_scattering.gaussian((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  assert g.n_terms() == 5
  assert approx_equal(g.array_of_a(),(1,-2,3,-4,5))
  assert approx_equal(g.array_of_b(),(-.1,.2,-.3,.4,-.5))
  assert approx_equal(g.c(), 6)
  assert approx_equal(g.at_stol_sq(3/4.), 13.4251206)
  assert approx_equal(g.at_stol(math.sqrt(3/4.)), 13.4251206)
  assert approx_equal(g.at_d_star_sq(3), 13.4251206)
  assert approx_equal(g.at_d_star(math.sqrt(3)), 13.4251206)
  s = pickle.dumps(g)
  l = pickle.loads(s)
  assert l.n_terms() == g.n_terms()
  assert approx_equal(l.array_of_a(), g.array_of_a())
  assert approx_equal(l.array_of_b(), g.array_of_b())
  assert approx_equal(l.c(), g.c())
  assert l.use_c()
  g = xray_scattering.gaussian((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5))
  s = pickle.dumps(g)
  l = pickle.loads(s)
  assert l.n_terms() == g.n_terms()
  assert approx_equal(l.array_of_a(), g.array_of_a())
  assert approx_equal(l.array_of_b(), g.array_of_b())
  assert approx_equal(l.c(), g.c())
  assert not l.use_c()

def exercise_it1992():
  e = xray_scattering.it1992("c1")
  assert e.table() == "IT1992"
  assert e.label() == "C"
  g = e.fetch()
  assert g.n_terms() == 4
  assert g.n_parameters() == 9
  assert approx_equal(g.array_of_a(), (2.31000, 1.02000, 1.58860, 0.865000))
  assert approx_equal(g.array_of_b(), (20.8439, 10.2075, 0.568700, 51.6512))
  assert approx_equal(g.c(), 0.215600)
  assert approx_equal(g.at_stol_sq(0), 5.99919997156)
  assert approx_equal(g.at_stol_sq(1./9), 2.26575563201)
  assert approx_equal(g.at_stol(1./9), 4.93537567523)
  assert approx_equal(g.at_d_star_sq(1./9), 4.04815863088)
  e = xray_scattering.it1992("yb2+", 1)
  assert e.label() == "Yb2+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 28.1209)
  assert approx_equal(g.array_of_b()[3], 20.3900)
  assert approx_equal(g.c(), 3.70983)
  e = xray_scattering.it1992("  YB3+")
  assert e.label() == "Yb3+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 27.8917)
  n = 0
  for e in xray_scattering.it1992_iterator():
    n += 1
    if (n == 215):
      assert e.label() == "Cf"
    d = xray_scattering.it1992(e.label(), 1)
    assert d.label() == e.label()
  assert n == 215
  i = xray_scattering.it1992_iterator()
  j = iter(i)
  assert i is j

def exercise_wk1995():
  e = xray_scattering.wk1995("c1")
  assert e.table() == "WK1995"
  assert e.label() == "C"
  g = e.fetch()
  assert approx_equal(g.array_of_a(),
    (2.657506,1.078079,1.490909,-4.241070,0.713791))
  assert approx_equal(g.array_of_b(),
    (14.780758,0.776775,42.086842,-0.000294,0.239535))
  assert approx_equal(g.c(), 4.297983)
  assert approx_equal(g.at_stol_sq(0), 5.99719834328)
  assert approx_equal(g.at_stol_sq(1./9), 2.26895371584)
  assert approx_equal(g.at_stol(1./9), 4.93735084739)
  assert approx_equal(g.at_d_star_sq(1./9), 4.04679561237)
  e = xray_scattering.wk1995("yb2+", 1)
  assert e.label() == "Yb2+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 28.443794)
  assert approx_equal(g.array_of_b()[4], 0.001463)
  assert approx_equal(g.c(), -23.214935)
  e = xray_scattering.wk1995("  YB3+")
  assert e.label() == "Yb3+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 28.191629)
  n = 0
  for e in xray_scattering.wk1995_iterator():
    n += 1
    if (n == 215):
      assert e.label() == "Pu6+"
    d = xray_scattering.wk1995(e.label(), 1)
    assert d.label() == e.label()
  assert n == 215
  i = xray_scattering.wk1995_iterator()
  j = iter(i)
  assert i is j

def ensure_common_symbols():
  lbl_it = []
  for e in xray_scattering.it1992_iterator(): lbl_it.append(e.label())
  lbl_it.sort()
  lbl_wk = []
  for e in xray_scattering.wk1995_iterator(): lbl_wk.append(e.label())
  lbl_wk.sort()
  assert lbl_it == lbl_wk

def ensure_correct_element_symbol():
  from cctbx.eltbx import tiny_pse
  for e in xray_scattering.it1992_iterator():
    l = e.label()
    if (l == "Cval"):
      s = "C"
    else:
      s = l[:2]
      if (len(s) > 1 and s[1] not in string.letters):
        s = s[:1]
    assert tiny_pse.table(l).symbol() == s
    assert tiny_pse.table(l.lower()).symbol() == s
    assert tiny_pse.table(l.upper()).symbol() == s

def run():
  exercise_gaussian()
  exercise_it1992()
  exercise_wk1995()
  ensure_common_symbols()
  ensure_correct_element_symbol()
  print "OK"

if (__name__ == "__main__"):
  run()
