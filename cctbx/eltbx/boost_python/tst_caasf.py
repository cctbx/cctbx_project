from cctbx.eltbx import caasf
from scitbx.test_utils import approx_equal

def exercise_it1992():
  c = caasf.it1992("const")
  assert c.table() == "IT1992"
  assert c.label() == "const"
  assert c.a() == (0,0,0,0)
  assert c.b() == (0,0,0,0)
  assert approx_equal(c.c(), 1)
  assert approx_equal(c.at_stol_sq(0), 1)
  assert approx_equal(c.at_stol_sq(1./9), 1)
  assert approx_equal(c.at_stol(1./9), 1)
  assert approx_equal(c.at_d_star_sq(1./9), 1)
  c = caasf.it1992("custom")
  assert c.table() == "IT1992"
  assert c.label() == "custom"
  assert c.a() == (0,0,0,0)
  assert c.b() == (0,0,0,0)
  assert approx_equal(c.c(), 0)
  c = caasf.it1992("c1")
  assert c.label() == "C"
  assert approx_equal(c.a(), (2.31000, 1.02000, 1.58860, 0.865000))
  assert approx_equal(c.b(), (20.8439, 10.2075, 0.568700, 51.6512))
  assert approx_equal(c.c(), 0.215600)
  c = caasf.it1992("yb2+", 1)
  assert c.label() == "Yb2+"
  assert approx_equal(c.a()[0], 28.1209)
  assert approx_equal(c.b()[3], 20.3900)
  assert approx_equal(c.c(), 3.70983)
  c = caasf.it1992("  YB3+")
  assert c.label() == "Yb3+"
  assert approx_equal(c.a()[0], 27.8917)
  n = 0
  for c in caasf.it1992_iterator():
    n += 1
    if (n == 1):
      assert c.label() == "const"
    elif (n == 2):
      assert c.label() == "custom"
    elif (n == 217):
      assert c.label() == "Cf"
    d = caasf.it1992(c.label(), 1)
    assert d.label() == c.label()
  assert n == 217
  i = caasf.it1992_iterator()
  j = iter(i)
  assert i is j

def exercise_wk1995():
  c = caasf.wk1995("const")
  assert c.table() == "WK1995"
  assert c.label() == "const"
  assert c.a() == (0,0,0,0,0)
  assert c.b() == (0,0,0,0,0)
  assert approx_equal(c.c(), 1)
  assert approx_equal(c.at_stol_sq(0), 1)
  assert approx_equal(c.at_stol_sq(1./9), 1)
  assert approx_equal(c.at_stol(1./9), 1)
  assert approx_equal(c.at_d_star_sq(1./9), 1)
  c = caasf.wk1995("custom")
  assert c.table() == "WK1995"
  assert c.label() == "custom"
  assert c.a() == (0,0,0,0,0)
  assert c.b() == (0,0,0,0,0)
  assert approx_equal(c.c(), 0)
  c = caasf.wk1995("c1")
  assert c.label() == "C"
  assert approx_equal(c.a(), (2.657506,1.078079,1.490909,-4.241070,0.713791))
  assert approx_equal(c.b(), (14.780758,0.776775,42.086842,-0.000294,0.239535))
  assert approx_equal(c.c(), 4.297983)
  c = caasf.wk1995("yb2+", 1)
  assert c.label() == "Yb2+"
  assert approx_equal(c.a()[0], 28.443794)
  assert approx_equal(c.b()[4], 0.001463)
  assert approx_equal(c.c(), -23.214935)
  c = caasf.wk1995("  YB3+")
  assert c.label() == "Yb3+"
  assert approx_equal(c.a()[0], 28.191629)
  n = 0
  for c in caasf.wk1995_iterator():
    n += 1
    if (n == 1):
      assert c.label() == "const"
    elif (n == 2):
      assert c.label() == "custom"
    elif (n == 217):
      assert c.label() == "Pu6+"
    d = caasf.wk1995(c.label(), 1)
    assert d.label() == c.label()
  assert n == 217, n
  i = caasf.wk1995_iterator()
  j = iter(i)
  assert i is j

def ensure_common_symbols():
  lbl_it = []
  for c in caasf.it1992_iterator(): lbl_it.append(c.label())
  lbl_it.sort()
  lbl_wk = []
  for c in caasf.wk1995_iterator(): lbl_wk.append(c.label())
  lbl_wk.sort()
  assert lbl_it == lbl_wk

def exercise_custom():
  c = caasf.custom("dummy", (1,-2,3,-4,5), (-.1,.2,-.3,.4,-.5), 6)
  assert c.label() == "dummy"
  assert c.n_ab() == 5
  assert approx_equal(c.a(),(1,-2,3,-4,5))
  assert approx_equal(c.b(),(-.1,.2,-.3,.4,-.5))
  assert approx_equal(c.c(), 6)
  c = caasf.custom("d", (), (), -2)
  assert c.label() == "d"
  assert c.n_ab() == 0
  assert c.a() == ()
  assert c.b() == ()
  assert approx_equal(c.c(), -2)

def run():
  exercise_it1992()
  exercise_wk1995()
  ensure_common_symbols()
  exercise_custom()
  print "OK"

if (__name__ == "__main__"):
  run()
