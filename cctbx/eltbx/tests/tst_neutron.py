from cctbx.eltbx import neutron
from libtbx.test_utils import approx_equal

def exercise():
  t = neutron.neutron_news_1992_table("eu")
  assert t.label() == "Eu"
  l = t.bound_coh_scatt_length()
  assert approx_equal(l.real, 7.22)
  assert approx_equal(l.imag, -1.26)
  assert approx_equal(t.abs_cross_sect(), 4530.)
  n = 0
  for t in neutron.neutron_news_1992_table_iterator():
    n += 1
    if (n == 1):
      assert t.label() == "H"
    elif (n == 86):
      assert t.label() == "U"
    u = neutron.neutron_news_1992_table(t.label())
    assert u.label() == t.label()
  assert n == 86

def run():
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run()
