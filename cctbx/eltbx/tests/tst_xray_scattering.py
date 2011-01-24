from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal
import pickle
import math

def exercise_basic():
  std_labels = xray_scattering.standard_labels_list()
  assert len(std_labels) == 217
  assert std_labels[:5] == ["H", "D", "T", "Hiso", "He"]
  assert std_labels[-1] == "Pu6+"
  for l in std_labels:
    assert xray_scattering.get_standard_label(
      label=l, exact=True, optional=False) == l
  assert xray_scattering.get_standard_label(label="na+") == "Na1+"
  assert xray_scattering.get_standard_label(label="na+") == "Na1+"
  assert xray_scattering.get_standard_label(label="o-") == "O1-"
  assert xray_scattering.get_standard_label(label="SI4+A") == "Si4+"
  assert xray_scattering.get_standard_label(label="SI1+") == "Si"
  assert xray_scattering.get_standard_label(label="SI1+",
    exact=True, optional=True) is None
  try:
    xray_scattering.get_standard_label(label="SI1+",
      exact=True, optional=False)
  except RuntimeError, e:
    assert str(e) == 'Unknown scattering type label: "SI1+"'
  else: raise Exception_expected
  #
  from cctbx.eltbx import tiny_pse
  for sl in std_labels:
    e, c = xray_scattering.get_element_and_charge_symbols(scattering_type=sl)
    assert e == "T" or tiny_pse.table(e, True).symbol() == e
    if (c != ""):
      assert len(c) == 2
      assert "123456789".find(c[0]) >= 0
      assert c[1] in ["+", "-"]

def exercise_gaussian():
  g = xray_scattering.gaussian(0)
  assert g.n_terms() == 0
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  assert g.n_parameters() == 1
  g = xray_scattering.gaussian(0, False)
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
  g = xray_scattering.gaussian(flex.double((1,2,3,4)))
  assert approx_equal(g.array_of_a(), (1,3))
  assert approx_equal(g.array_of_b(), (2,4))
  assert approx_equal(g.c(), 0)
  assert not g.use_c()
  g = xray_scattering.gaussian(flex.double((1,2,3,4)), 0, True)
  assert approx_equal(g.c(), 0)
  assert g.use_c()
  g = xray_scattering.gaussian(flex.double((1,2,3,4)), 5)
  assert approx_equal(g.c(), 5)
  assert g.use_c()
  g = xray_scattering.gaussian((1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  assert g.n_terms() == 5
  assert approx_equal(g.array_of_a(),(1,-2,3,-4,5))
  assert approx_equal(g.array_of_b(),(-.1,.2,-.3,.4,-.5))
  assert approx_equal(g.c(), 6)
  assert approx_equal(g.at_stol_sq(3/4.), 13.4251206)
  assert approx_equal(g.at_stol(math.sqrt(3/4.)), 13.4251206)
  assert approx_equal(g.at_d_star_sq(3), 13.4251206)
  assert approx_equal(g.at_d_star_sq(flex.double([3,4])),
    [13.4251206, 15.079612])
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
  #
  g = xray_scattering.gaussian(*zip(*[
    (2.51340127252, 31.8053433708),
    (1.74867019409, 0.445605499982),
    (1.72398202356, 10.5831679451)])) # C 3-gaussian
  assert not g.use_c()
  e = g.electron_density
  assert approx_equal(e(r=0, b_iso=0), 264.731533932)
  assert approx_equal(e(r=0, b_iso=1), 47.3615000971)
  assert approx_equal(e(r=1, b_iso=0), 0.233911429529)
  assert approx_equal(e(r=1, b_iso=1), 0.243343387016)
  g = xray_scattering.gaussian(
    (2.31000, 1.02000, 1.58860, 0.865000),
    (20.8439, 10.2075, 0.568700, 51.6512),
    0.215600) # C it1992
  assert g.use_c()
  e = g.electron_density
  assert approx_equal(e(r=0, b_iso=0.1), 435.677592698)
  assert approx_equal(e(r=0, b_iso=1), 47.9420591405)
  assert approx_equal(e(r=1, b_iso=0.1), 0.241082494004)
  assert approx_equal(e(r=1, b_iso=1), 0.248806720643)

def exercise_n_gaussian():
  assert xray_scattering.n_gaussian_table_size() == 213
  assert xray_scattering.n_gaussian_table_index("H") == 0
  assert xray_scattering.n_gaussian_table_index("Pu6+") == 212
  for n_terms in [6,5,4,3,2,1]:
    e = xray_scattering.n_gaussian_table_entry(0, n_terms)
    assert e.label() == "H"
    g = e.gaussian()
    assert g.n_terms() == n_terms
    assert approx_equal(g.at_x(0), 1, eps=0.01+1.e-6)
    assert e.max_stol() > 0
    assert e.d_min() > 0
    assert e.max_relative_error() > 0
  for i_entry in xrange(xray_scattering.n_gaussian_table_size()):
    for n_terms in [6,5,4,3,2,1]:
      e = xray_scattering.n_gaussian_table_entry(i_entry, n_terms)
      assert e.gaussian().n_terms() == n_terms
      f = xray_scattering.n_gaussian_table_entry(e.label(), n_terms)
      assert f.label() == e.label()
      assert f.gaussian().n_terms() == n_terms
    for d_min in [0,10]:
      for max_relative_error in [0,0.5]:
        e = xray_scattering.n_gaussian_table_entry(
          i_entry, d_min, max_relative_error)
        f = xray_scattering.n_gaussian_table_entry(
          e.label(), d_min, max_relative_error)
        assert f.label() == e.label()
        if (d_min == 0):
          n_terms = 6
        else:
          n_terms = 1
        assert e.gaussian().n_terms() == n_terms
        assert f.gaussian().n_terms() == n_terms
  label = "Be"
  be_max_stols = []
  be_max_relative_errors = []
  for n_terms in [6,5,4,3,2,1]:
    e = xray_scattering.n_gaussian_table_entry(label, n_terms)
    be_max_stols.append(e.max_stol())
    be_max_relative_errors.append(e.max_relative_error())
  for n_terms,stol,max_relative_error in zip([6,5,4,3,2,1],
                                             be_max_stols,
                                             be_max_relative_errors):
    e = xray_scattering.n_gaussian_table_entry(
      "Be", 1/(2*stol)+1.e-6, max_relative_error+1.e-6)
    assert e.gaussian().n_terms() == n_terms
    assert approx_equal(e.max_stol(), stol)
    if (n_terms < 6):
      e = xray_scattering.n_gaussian_table_entry(
        "Be", 1/(2*stol)-1.e-6, max(be_max_relative_errors)+1.e-6)
      assert e.gaussian().n_terms() == min(n_terms+1, 6)
  assert be_max_relative_errors[1] > be_max_relative_errors[2]

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
  e = xray_scattering.it1992("yb2+", True)
  assert e.label() == "Yb2+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 28.1209)
  assert approx_equal(g.array_of_b()[3], 20.3900)
  assert approx_equal(g.c(), 3.70983)
  e = xray_scattering.it1992("  yB3+")
  assert e.label() == "Yb3+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 27.8917)
  n = 0
  for e in xray_scattering.it1992_iterator():
    n += 1
    if (n == 213):
      assert e.label() == "Cf"
    else:
      assert e.label() != "Cf"
    d = xray_scattering.it1992(e.label(), True)
    assert d.label() == e.label()
  assert n == 213
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
  e = xray_scattering.wk1995("yb2+", True)
  assert e.label() == "Yb2+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 28.443794)
  assert approx_equal(g.array_of_b()[4], 0.001463)
  assert approx_equal(g.c(), -23.214935)
  e = xray_scattering.wk1995("  yB3+")
  assert e.label() == "Yb3+"
  g = e.fetch()
  assert approx_equal(g.array_of_a()[0], 28.191629)
  n = 0
  for e in xray_scattering.wk1995_iterator():
    n += 1
    if (n == 213):
      assert e.label() == "Pu6+"
    else:
      assert e.label() != "Pu6+"
    d = xray_scattering.wk1995(e.label(), True)
    assert d.label() == e.label()
  assert n == 213
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
  assert lbl_wk == lbl_it
  lbl_ng = []
  for i_entry in xrange(xray_scattering.n_gaussian_table_size()):
    lbl_ng.append(xray_scattering.n_gaussian_table_entry(i_entry, 6).label())
  lbl_ng.sort()
  assert lbl_ng == lbl_it
  #
  for label in xray_scattering.standard_labels_list():
    it = xray_scattering.it1992(label, True).fetch()
    wk = xray_scattering.wk1995(label, True).fetch()
    ng = xray_scattering.n_gaussian_table_entry(label, 0, 0).gaussian()
    assert approx_equal(wk.at_stol(0)/it.at_stol(0), 1, 5.e-3)

def ensure_correct_element_symbol():
  from cctbx.eltbx import tiny_pse
  for e in xray_scattering.it1992_iterator():
    l = e.label()
    e, c = xray_scattering.get_element_and_charge_symbols(
      scattering_type=l, exact=False)
    assert tiny_pse.table(l).symbol() == e
    assert tiny_pse.table(l.lower()).symbol() == e
    assert tiny_pse.table(l.upper()).symbol() == e

def run():
  exercise_basic()
  exercise_gaussian()
  exercise_n_gaussian()
  exercise_it1992()
  exercise_wk1995()
  ensure_common_symbols()
  ensure_correct_element_symbol()
  print "OK"

if (__name__ == "__main__"):
  run()
