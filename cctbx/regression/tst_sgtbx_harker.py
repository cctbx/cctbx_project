from __future__ import absolute_import, division, print_function
from cctbx.sgtbx import harker
from cctbx import crystal
from cctbx import sgtbx
from libtbx.test_utils import approx_equal as ae

def run():
  uc = (10,12,14,90,90,90)
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P 1")
  assert harker.planes_cartesian(cs).count() == 0
  assert harker.planes_cartesian(cs).min_distance((0,0,0)) is None
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P 2 2 2")
  assert harker.planes_cartesian(cs).count() == 3
  assert ae(harker.planes_cartesian(cs).min_distance((0.123,0.234,0)), 0)
  assert ae(harker.planes_cartesian(cs).min_distance((0.123,0,0.234)), 0)
  assert ae(harker.planes_cartesian(cs).min_distance((0,0.123,0.234)), 0)
  assert ae(harker.planes_cartesian(cs).min_distance((0.1,0.234,0.345)), 1)
  assert ae(harker.planes_cartesian(cs).min_distance((0.234,0.1,0.345)), 1.2)
  assert ae(harker.planes_cartesian(cs).min_distance((0.234,0.345,0.1)), 1.4)
  assert ae(harker.planes_cartesian(cs).min_distance((0.2,0.2,0.1)), 1.4)
  assert ae(harker.planes_cartesian(cs).min_distance((0.2,0.1,0.1)), 1.2)
  assert ae(harker.planes_cartesian(cs).min_distance((0.1,0.2,0.2)), 1)
  assert ae(harker.planes_cartesian(cs).min_distance((-0.2,0.1,0.1)), 1.2)
  assert ae(harker.planes_cartesian(cs).min_distance((0.2,-0.1,0.1)), 1.2)
  assert ae(harker.planes_cartesian(cs).min_distance((0.2,0.1,-0.1)), 1.2)
  assert ae(harker.planes_cartesian(cs).min_distance((-0.2,-0.1,-0.1)), 1.2)
  uc = (10,10,12,90,90,120)
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P 3")
  assert harker.planes_cartesian(cs).count() == 1
  assert ae(harker.planes_cartesian(cs).min_distance((0,0,0)), 0)
  assert ae(harker.planes_cartesian(cs).min_distance((0.1,0.2,1/3.)), 4)
  assert ae(harker.planes_cartesian(cs).min_distance((0.1,0.2,1/2.)), 6)
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P 31")
  assert harker.planes_cartesian(cs).count() == 1
  assert ae(harker.planes_cartesian(cs).min_distance((0,0,0)), 4)
  assert ae(harker.planes_cartesian(cs).min_distance((0.1,0.2,1/3.)), 0)
  assert ae(harker.planes_cartesian(cs).min_distance((0.1,0.2,1/2.)), 2)
  uc = (10,10,12,90,90,90)
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P 4 2 2")
  assert harker.planes_cartesian(cs).count() == 5
  cs = crystal.symmetry(unit_cell=uc, space_group_symbol="P 41 2 2")
  assert harker.planes_cartesian(cs).count() == 6
  ps = crystal.symmetry(unit_cell=uc, space_group_symbol="P 4/m m m")
  for x,d,t in (((0.1,0.2,0.3),0.6,"1"),
                ((0.1,0.2,0.75),0,"1"),
                ((0.1,0.0,0.3),0,"m")):
    ss = crystal.special_position_settings(ps).site_symmetry(x)
    assert ss.point_group_type() == t
    eq = sgtbx.sym_equiv_sites(ss)
    for x in eq.coordinates():
      assert ae(harker.planes_cartesian(cs).min_distance(x), d)
  print("OK")

if (__name__ == "__main__"):
  run()
