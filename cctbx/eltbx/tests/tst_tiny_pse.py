from __future__ import absolute_import, division, print_function
from cctbx.eltbx import tiny_pse
from libtbx.test_utils import approx_equal

def exercise():
  t = tiny_pse.table("SI")
  assert t.atomic_number() == 14
  assert t.symbol() == "Si"
  assert t.name() == "silicon"
  assert approx_equal(t.weight(), 28.086)
  n = 0
  for t in tiny_pse.table_iterator():
    n += 1
    if (n == 1):
      assert t.symbol() == "H"
    elif (n == 104):
      assert t.atomic_number() == 103
      assert t.symbol() == "Lr"
    u = tiny_pse.table(t.symbol())
    assert u.symbol() == t.symbol()
  assert n == 104

def run():
  exercise()
  print("OK")

if (__name__ == "__main__"):
  run()
