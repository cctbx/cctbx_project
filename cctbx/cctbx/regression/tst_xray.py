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
  ys = xs.deep_copy_scatterers()
  ys.add_scatterers(ys.scatterers())
  assert ys.scatterers().size() == 4
  assert tuple(ys.special_position_indices()) == (0, 1, 2, 3)
  ys.add_scatterer(ys.scatterers()[0])
  assert ys.scatterers().size() == 5
  assert tuple(ys.special_position_indices()) == (0, 1, 2, 3, 4)
  sx = xs.primitive_setting()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 62 2 2"
  sx = xs.change_hand()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 64 2 2"
  assert approx_equal(sx.scatterers()[0].site, (-1./2, -1./2, -1./3))
  assert approx_equal(sx.scatterers()[1].site, (-0.19700, 0.19700, -0.833333))
  p1 = xs.expand_to_p1()
  assert p1.scatterers().size() == 9
  sh = p1.apply_shift((0.2,0.3,-1/6.))
  assert approx_equal(sh.scatterers()[0].site, (0.7,0.8,1/6.))
  assert approx_equal(sh.scatterers()[3].site, (0.3970,0.1030,2/3.))
  sl = sh[:1]
  assert sl.scatterers().size() == 1
  assert sl.scatterers()[0].label == sh.scatterers()[0].label
  sl = sh[1:4]
  assert sl.scatterers().size() == 3
  for i in xrange(3):
    assert sl.scatterers()[i].label == sh.scatterers()[i+1].label
  xs.scatterers().set_occupancies(flex.double((0.5,0.2)))
  s = xs.sort(by_value="occupancy")
  assert approx_equal(s.scatterers().extract_occupancies(), (0.2,0.5))
  assert s.scatterers()[0].label == "O1"
  assert s.scatterers()[1].label == "Si1"

def run():
  exercise_structure()
  print "OK"

if (__name__ == "__main__"):
  run()
