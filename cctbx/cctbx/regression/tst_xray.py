from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_structure():
  cs = crystal.symmetry((5.01, 5.01, 5.47, 90, 90, 120), "P 62 2 2")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", (1./2, 1./2, 1./3)),
    xray.scatterer("O1", (0.19700, -0.19700, 0.83333))))
  xs = xray.structure(sp, scatterers)
  assert xs.scatterers().size() == 2
  assert tuple(xs.special_position_indices()) == (0, 1)
  xs.all_apply_symmetry()
  assert tuple(xs.special_position_indices()) == (0, 1)
  sx = xs.change_hand()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 64 2 2"
  assert approx_equal(sx.scatterers()[0].site, (-1./2, -1./2, -1./3))
  assert approx_equal(sx.scatterers()[1].site, (-0.19700, 0.19700, -0.833333))

def run():
  exercise_structure()
  print "OK"

if (__name__ == "__main__"):
  run()
