from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.development import debug_utils
from libtbx.test_utils import approx_equal
import sys

def exercise_symmetry():
  xs = crystal.symmetry()
  xs = crystal.symmetry(uctbx.unit_cell((3,4,5)))
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  xs = crystal.symmetry("3,4,5", "P 2 2 2")
  xs = crystal.symmetry([4,5,6], space_group_info=xs.space_group_info())
  xs = crystal.symmetry([4,5,6], space_group=xs.space_group())
  xs = crystal.symmetry([4,5,6], space_group=str(xs.space_group_info()))
  assert str(xs.unit_cell()) == "(4, 5, 6, 90, 90, 90)"
  assert xs.unit_cell().is_similar_to(uctbx.unit_cell((4,5,6)))
  assert str(xs.space_group_info()) == "P 2 2 2"
  assert xs.space_group() == xs.space_group_info().group()
  assert xs.is_compatible_unit_cell()
  try: xs = crystal.symmetry((3,4,5), "P 4 2 2")
  except: pass
  else: raise AssertionError, "Exception expected."
  xs = crystal.symmetry(
    (3,4,5), "P 4 2 2", assert_is_compatible_unit_cell=False,
    force_compatible_unit_cell=False)
  assert not xs.is_compatible_unit_cell()
  xs = crystal.symmetry(
    (3,4,5), "P 4 2 2", assert_is_compatible_unit_cell=False)
  assert xs.is_compatible_unit_cell()
  xs = crystal.symmetry((5,5,29,90,90,120), "R 3")
  ps = xs.primitive_setting()
  assert ps.unit_cell().is_similar_to(
    uctbx.unit_cell((10.0885, 10.0885, 10.0885, 28.6956, 28.6956, 28.6956)))
  assert str(ps.space_group_info()) == "R 3 :R"
  rs = ps.as_reference_setting()
  assert rs.unit_cell().is_similar_to(xs.unit_cell())
  assert str(rs.space_group_info()) == "R 3 :H"
  cb = xs.change_of_basis_op_to_minimum_cell()
  mc = xs.minimum_cell()
  cm = mc.change_basis(cb.inverse())
  assert cm.unit_cell().is_similar_to(xs.unit_cell())
  assert cm.space_group() == xs.space_group()
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
  xs = crystal.symmetry((5,3,4), "P 2 2 2")
  p1 = xs.cell_equivalent_p1()
  assert p1.unit_cell().is_similar_to(uctbx.unit_cell((5,3,4)))
  assert p1.space_group().order_z() == 1
  ps = xs.patterson_symmetry()
  assert ps.unit_cell().is_similar_to(xs.unit_cell())
  assert str(ps.space_group_info()) == "P m m m"
  bc = ps.best_cell()
  assert bc.unit_cell().is_similar_to(uctbx.unit_cell((3,4,5)))
  assert str(bc.space_group_info()) == "P m m m"
  xs = crystal.symmetry((5,3,4,90,130,90), "P 1 2 1")
  bc = xs.best_cell()
  assert bc.unit_cell().is_similar_to(
    uctbx.unit_cell((3.91005,3,4,90,101.598,90)))
  assert str(bc.space_group_info()) == "P 1 2 1"
  cb = xs.change_of_basis_op_to_best_cell()
  assert str(cb.c()) in ["x,-y,x-z", "-x,-y,-x+z"]
  assert bc.change_basis("x,-y,x-z").unit_cell().is_similar_to(
         bc.change_basis("-x,-y,-x+z").unit_cell())
  asu = xs.direct_space_asu()
  assert asu.hall_symbol == " P 2y"
  assert len(asu.facets) == 6
  assert asu.unit_cell is xs.unit_cell()
  asu_mappings = xs.asu_mappings(buffer_thickness=2.364)
  assert approx_equal(asu_mappings.buffer_thickness(), 2.364)

def exercise_special_position_settings():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  sp = crystal.special_position_settings(xs, 1, 2, True, False)
  assert sp.min_distance_sym_equiv() == 1
  assert sp.u_star_tolerance() == 2
  assert sp.assert_is_positive_definite() == True
  assert sp.assert_min_distance_sym_equiv() == False
  assert sp.site_symmetry((0,0,0)).multiplicity() == 1
  assert str(sp.sym_equiv_sites((0,0,0)).special_op()) == "0,0,0"

def exercise_site_symmetry(space_group_info):
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry=crystal.symmetry(
      unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
      space_group_info=space_group_info))
  z2p_op = space_group_info.group().z2p_op()
  special_position_settings_p = crystal.special_position_settings(
    crystal_symmetry=special_position_settings.change_basis(z2p_op),
    min_distance_sym_equiv
      =special_position_settings.min_distance_sym_equiv()*0.99)
  wyckoff_table = space_group_info.wyckoff_table()
  for i_position in xrange(wyckoff_table.size()):
    site_symmetry = wyckoff_table.random_site_symmetry(
      special_position_settings=special_position_settings,
      i_position=i_position)
    s = site_symmetry.special_op()
    assert s.multiply(s) == s
    for m in site_symmetry.matrices():
      assert m.multiply(s) == s
    tab = sgtbx.site_symmetry_table()
    tab.process(site_symmetry)
    ss_ops = tab.get(0)
    assert ss_ops.multiplicity() == site_symmetry.multiplicity()
    assert ss_ops.multiplicity() * ss_ops.n_matrices() \
        == site_symmetry.space_group().order_z()
    site_p = z2p_op.c() * site_symmetry.exact_site()
    site_symmetry_p = special_position_settings_p.site_symmetry(site_p)
    ss_ops_p = ss_ops.change_basis(z2p_op)
    assert ss_ops_p.multiplicity() == site_symmetry_p.multiplicity()
    assert ss_ops_p.special_op() == site_symmetry_p.special_op()
    assert ss_ops_p.multiplicity() * ss_ops_p.n_matrices() \
        == site_symmetry_p.space_group().order_z()
    references = [str(m) for m in site_symmetry_p.matrices()]
    testees = [str(m) for m in ss_ops_p.matrices()]
    references.sort()
    testees.sort()
    assert " ".join(testees) == " ".join(references)

def run_call_back(flags, space_group_info):
  exercise_site_symmetry(space_group_info)

def run():
  exercise_symmetry()
  exercise_special_position_settings()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
