from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex

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

def run():
  exercise_structure()
  print "OK"

if (__name__ == "__main__"):
  run()
