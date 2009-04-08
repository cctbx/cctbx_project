from cctbx.eltbx import covalent_radii
from libtbx.test_utils import approx_equal

def exercise():
  t = covalent_radii.table("sI")
  assert t.label() == "Si"
  assert approx_equal(t.radius(), 1.11)
  assert approx_equal(t.esd(), 0.02)
  n = 0
  for t in covalent_radii.table_iterator():
    n += 1
    if (n == 1):
      assert t.label() == "H"
    elif (n == 96):
      assert t.label() == "Cm"
    u = covalent_radii.table(t.label())
    assert u.label() == t.label()
  assert n == 96

def run():
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run()
