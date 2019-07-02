from __future__ import absolute_import, division, print_function
from cctbx.eltbx import sasaki
from libtbx.test_utils import approx_equal

def verify(table, wave_length, fp, fdp):
  fp_fdp = table.at_angstrom(wave_length)
  assert approx_equal(fp_fdp.fp(), fp)
  assert approx_equal(fp_fdp.fdp(), fdp)

def exercise():
  t = sasaki.table("SI")
  assert t.label() == "Si"
  assert t.atomic_number() == 14
  f = t.at_angstrom(2)
  assert f.is_valid_fp()
  assert f.is_valid_fdp()
  assert f.is_valid()
  from cctbx import factor_kev_angstrom
  assert approx_equal(f.fp(), t.at_kev(factor_kev_angstrom / 2).fp())
  assert approx_equal(f.fdp(), t.at_kev(factor_kev_angstrom / 2).fdp())
  assert approx_equal(f.fp(), t.at_ev(1000 * factor_kev_angstrom / 2).fp())
  assert approx_equal(f.fdp(), t.at_ev(1000 * factor_kev_angstrom / 2).fdp())
  c = f.as_complex()
  assert c.real == f.fp()
  assert c.imag == f.fdp()
  verify(t, 0.1, -0.0226, 0.0010) # wide first
  verify(t, 2.89, 0.3824, 1.0517) # wide last
  verify(t, 6.7289, -6.9495, 4.1042) # K edge
  t = sasaki.table("ag")
  assert t.label() == "Ag"
  assert t.atomic_number() == 47
  verify(t, 3.2426, -8.6870, 14.1062) # L1 edge
  verify(t, 3.5169, -16.5123, 13.7723) # L2 edge
  verify(t, 3.6995, -29.9836, 10.5137) # L3 edge, right before edge
  verify(t, 3.6996, -32.2967, 3.1759) # L3 edge, right after edge
  n = 0
  for t in sasaki.table_iterator():
    n += 1
    if (n == 1):
      assert t.label() == "Li"
    elif (n == 82):
      assert t.label() == "U"
    u = sasaki.table(t.label())
    assert u.label() == t.label()
  assert n == 82

def run():
  exercise()
  print("OK")

if (__name__ == "__main__"):
  run()
