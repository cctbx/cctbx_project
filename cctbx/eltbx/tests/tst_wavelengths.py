from __future__ import absolute_import, division, print_function
from cctbx.eltbx import wavelengths
from libtbx.test_utils import approx_equal

def exercise():
  from cctbx import factor_kev_angstrom
  w = wavelengths.characteristic("CU")
  assert w.label() == "Cu"
  assert approx_equal(w.as_angstrom(), 1.5418)
  assert approx_equal(w.as_kev(), factor_kev_angstrom / 1.5418)
  assert approx_equal(w.as_ev() / 1000, factor_kev_angstrom / 1.5418)
  n = 0
  for w in wavelengths.characteristic_iterator():
    n += 1
    uu = wavelengths.characteristic(w.label())
    assert uu.label() == w.label()
    assert uu.as_ev() == w.as_ev()
  assert n == 15

def run():
  exercise()
  print("OK")

if (__name__ == "__main__"):
  run()
