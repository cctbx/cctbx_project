from cctbx.eltbx import icsd_radii
from libtbx.test_utils import approx_equal

def exercise():
  t = icsd_radii.table("sI4+")
  assert t.label() == "Si4+"
  assert approx_equal(t.radius(), 0.26)
  n = 0
  for t in icsd_radii.table_iterator():
    n += 1
    if (n == 1):
      assert t.label() == "H"
    elif (n == 442):
      assert t.label() == "Lr"
    u = icsd_radii.table(t.label())
    assert u.label() == t.label()
  assert n == 442

def run():
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run()
