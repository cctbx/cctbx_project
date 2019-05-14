from __future__ import absolute_import, division, print_function
from cctbx.eltbx import henke
from libtbx.test_utils import approx_equal

def verify(table, energy, fp, fdp):
  fp_fdp = table.at_ev(energy)
  assert approx_equal(fp_fdp.fp(), fp)
  assert approx_equal(fp_fdp.fdp(), fdp)

def exercise():
  t = henke.table("SI")
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
  verify(t, 10.0, -9999.00, 4.00688)
  verify(t, 29.3, 4.04139-14, 0.371742)
  verify(t, 30000.0, 14.0266-14, 0.0228459)
  n = 0
  for t in henke.table_iterator():
    n += 1
    if (n == 1):
      assert t.label() == "H"
    elif (n == 92):
      assert t.label() == "U"
    u = henke.table(t.label())
    assert u.label() == t.label()
  assert n == 92

def run():
  exercise()
  print("OK")

if (__name__ == "__main__"):
  run()
