from cctbx import uctbx
from cctbx import crystal

def exercise_symmetry():
  xs = crystal.symmetry()
  xs = crystal.symmetry(uctbx.unit_cell((3,4,5)))
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  xs = crystal.symmetry("3,4,5", "P 2 2 2")
  xs = crystal.symmetry([4,5,6], space_group_info=xs.space_group_info())
  xs = crystal.symmetry([4,5,6], space_group=xs.space_group())
  assert str(xs.unit_cell()) == "(4, 5, 6, 90, 90, 90)"
  assert xs.unit_cell().is_similar_to(uctbx.unit_cell((4,5,6)))
  assert str(xs.space_group_info()) == "P 2 2 2"
  assert xs.space_group() == xs.space_group_info().group()
  assert xs.is_compatible_unit_cell()
  try: xs = crystal.symmetry((3,4,5), "P 4 2 2")
  except: pass
  else: raise AssertionError, "Exception expected."
  xs = crystal.symmetry(
    (3,4,5), "P 4 2 2", assert_is_compatible_unit_cell=00000,
    force_compatible_unit_cell=00000)
  assert not xs.is_compatible_unit_cell()
  xs = crystal.symmetry(
    (3,4,5), "P 4 2 2", assert_is_compatible_unit_cell=00000)
  assert xs.is_compatible_unit_cell()
  xs = crystal.symmetry((5,5,29,90,90,120), "R 3")
  ps = xs.primitive_setting()
  assert ps.unit_cell().is_similar_to(
    uctbx.unit_cell((10.0885, 10.0885, 10.0885, 28.6956, 28.6956, 28.6956)))
  assert str(ps.space_group_info()) == "R 3 :R"
  rs = ps.as_reference_setting()
  assert rs.unit_cell().is_similar_to(xs.unit_cell())
  assert str(rs.space_group_info()) == "R 3 :H"
  cb = xs.change_of_basis_op_to_niggli_cell()
  assert str(cb.c()) == "y-z,-x-z,3*z"
  nc = xs.niggli_cell()
  assert nc.unit_cell().is_similar_to(
    uctbx.unit_cell((5, 5, 10.0885, 75.6522, 75.6522, 60)))
  assert str(nc.space_group_info()) == "Hall:  R 3 (x+z,-y+z,-3*z)"
  assert nc.unit_cell().is_niggli_cell()
  cn = nc.change_basis(cb.inverse())
  assert cn.unit_cell().is_similar_to(xs.unit_cell())
  assert cn.space_group() == xs.space_group()
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  p1 = xs.cell_equivalent_p1()
  assert p1.unit_cell().is_similar_to(uctbx.unit_cell((3,4,5)))
  assert p1.space_group().order_z() == 1
  ps = xs.patterson_symmetry()
  assert ps.unit_cell().is_similar_to(xs.unit_cell())
  assert str(ps.space_group_info()) == "P m m m"

def exercise_special_position_settings():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  sp = crystal.special_position_settings(xs, 1, 2, 0001, 00000)
  assert sp.min_distance_sym_equiv() == 1
  assert sp.u_star_tolerance() == 2
  assert sp.assert_is_positive_definite() == 0001
  assert sp.assert_min_distance_sym_equiv() == 00000
  assert sp.site_symmetry((0,0,0)).multiplicity() == 1
  assert str(sp.sym_equiv_sites((0,0,0)).special_op()) == "0,0,0"

def run():
  exercise_symmetry()
  exercise_special_position_settings()
  print "OK"

if (__name__ == "__main__"):
  run()
