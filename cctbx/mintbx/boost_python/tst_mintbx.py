from cctbx import mintbx
from cctbx import uctbx
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_k_b_scaling():
  uc = uctbx.unit_cell((11,13,17))
  i = flex.miller_index(((1,2,3),(2,3,4)))
  m = flex.int((1,2))
  d = flex.double((1,2))
  kb = mintbx.k_b_scaling_target_and_gradients(
    uc, i, m, d, d, 1, 0, True, True)
  assert approx_equal(kb.target(), 0)
  assert approx_equal(kb.gradient_k(), 0)
  assert not kb.anisotropic_flag()
  assert approx_equal(kb.gradient_b_iso(), 0)
  kb = mintbx.k_b_scaling_target_and_gradients(
    uc, i, m, d, d, 1, (0,0,0,0,0,0), True, True)
  assert approx_equal(kb.target(), 0)
  assert approx_equal(kb.gradient_k(), 0)
  assert kb.anisotropic_flag()
  assert approx_equal(kb.gradients_b_cif(), (0,0,0,0,0,0))

def run():
  exercise_k_b_scaling()
  print "OK"

if (__name__ == "__main__"):
  run()
