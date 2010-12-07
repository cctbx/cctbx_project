from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from cStringIO import StringIO
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
  else: raise Exception_expected
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
  assert str(nc.space_group_info()) == "R 3 :H (x+z,-y+z,-3*z)"
  assert nc.unit_cell().is_niggli_cell()
  cn = nc.change_basis(cb.inverse())
  assert cn.unit_cell().is_similar_to(xs.unit_cell())
  assert cn.space_group() == xs.space_group()
  xs = crystal.symmetry((3,3,4,90,90,120), "P 31")
  ih = xs.inverse_hand()
  assert ih.unit_cell().is_similar_to(xs.unit_cell())
  assert str(ih.space_group_info()) == "P 32"
  xs = crystal.symmetry((5,3,4), "P 2 2 2")
  p1 = xs.cell_equivalent_p1()
  assert p1.unit_cell().is_similar_to(uctbx.unit_cell((5,3,4)))
  assert p1.space_group().order_z() == 1
  ri = xs.reflection_intensity_symmetry(anomalous_flag=True)
  assert ri.unit_cell().is_similar_to(xs.unit_cell())
  assert str(ri.space_group_info()) == "P 2 2 2"
  ri = xs.reflection_intensity_symmetry(anomalous_flag=False)
  assert ri.unit_cell().is_similar_to(xs.unit_cell())
  assert str(ri.space_group_info()) == "P m m m"
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
  assert len(asu.cuts) == 6
  assert asu.unit_cell is xs.unit_cell()
  asu_mappings = xs.asu_mappings(buffer_thickness=2.364)
  assert approx_equal(asu_mappings.buffer_thickness(), 2.364)
  assert approx_equal(xs.average_b_cart((1,2,3,4,5,6)), (1,2,3,0,5,0))
  #
  s = crystal.symmetry(
    unit_cell=(10,20,30,90,90,80),
    space_group_symbol="A 1 1 2")
  assert approx_equal(
    s.subtract_continuous_allowed_origin_shifts(translation_cart=[1,2,3]),
    [1,2,0])
  #
  for anomalous_flag in [False, True]:
    m = s.build_miller_set(
      anomalous_flag=anomalous_flag, d_min=8.1, d_max=8.5)
    assert m.anomalous_flag() == anomalous_flag
    if (not anomalous_flag):
      assert list(m.indices()) == [(1,0,2), (0,2,2)]
    else:
      assert list(m.indices()) == [(1,0,2), (-1,0,-2), (0,2,2), (0,-2,-2)]

def exercise_correct_rhombohedral_setting_if_necessary():
  for symbol in sgtbx.rhombohedral_hermann_mauguin_symbols:
    for p,z in [("20 20 21 90 90 120", "R"), ("31 31 31 85 85 85", "H")]:
      uc = uctbx.unit_cell(p)
      cs = crystal.symmetry(
        unit_cell=uc,
        space_group_symbol=symbol+":"+z,
        correct_rhombohedral_setting_if_necessary=True)
      assert cs.unit_cell().is_similar_to(uc)
      other_z = {
        "R": "H",
        "H": "R"}[z]
      assert not show_diff(
        cs.space_group_info().type().lookup_symbol(),
        symbol+" :"+other_z)
  cs = crystal.symmetry(
    unit_cell="20 20 21 90 89 120",
    space_group_symbol="R3:R",
    correct_rhombohedral_setting_if_necessary=True,
    assert_is_compatible_unit_cell=False)
  sio = StringIO()
  cs.show_summary(f=sio)
  assert not show_diff(sio.getvalue(), """\
Unit cell: (20.3388, 20.3388, 20.3388, 98.9315, 98.9315, 98.9315)
Space group: R 3 :R (No. 146)
""")
  cs = crystal.symmetry(
    unit_cell="31 31 31 85 85 86",
    space_group_symbol="R3:H",
    correct_rhombohedral_setting_if_necessary=True,
    assert_is_compatible_unit_cell=False)
  sio = StringIO()
  cs.show_summary(f=sio)
  assert not show_diff(sio.getvalue(), """\
Unit cell: (36.4146, 36.4146, 31, 90, 90, 120)
Space group: R 3 :H (No. 146)
""")

def exercise_select_crystal_symmetry():
  xs1 = crystal.symmetry(unit_cell   = "23,30,40,90,90,90",
                         space_group = "P212121" )
  xs2 = crystal.symmetry(unit_cell   = "20,30,40,90,90,90",
                         space_group = "P222" )
  resulting_symmetry = crystal.select_crystal_symmetry( from_command_line     = None,
                                                        from_parameter_file   = None,
                                                        from_coordinate_files = [xs1],
                                                        from_reflection_files = [xs2] )
  assert list( xs2.unit_cell().parameters()  ) == list( resulting_symmetry.unit_cell().parameters() )
  resulting_symmetry = crystal.select_crystal_symmetry( from_command_line     = None,
                                                        from_parameter_file   = None,
                                                        from_coordinate_files = [xs2],
                                                        from_reflection_files = [xs1] )
  assert list( xs1.unit_cell().parameters()  ) == list( resulting_symmetry.unit_cell().parameters() )

  resulting_symmetry = None
  try:
    resulting_symmetry = crystal.select_crystal_symmetry( from_command_line     = None,
                                                          from_parameter_file   = None,
                                                          from_coordinate_files = [None],
                                                          from_reflection_files = [None] )
  except AssertionError ,e :
    assert str(e)=="No unit cell and symmetry information supplied"
  else: raise Exception_expected

def verify_definitions_in_paper_zwart_2007():
  # Verification of definitions in Peter Zwart's paper for the
  # CCP4 Study Weekend Jan 2007.
  #
  cb_symbol_xyz = "x-y,x+y,z"
  cb_symbol_abc = "1/2*a-1/2*b,1/2*a+1/2*b,c"
  #
  # Verify the claim that cb_symbol_abc is the inverse transpose of
  # cb_symbol_xyz.
  cb_mx_xyz = sgtbx.rt_mx(cb_symbol_xyz, r_den=12, t_den=144)
  assert sgtbx.rt_mx(cb_mx_xyz.r().inverse().transpose()).as_xyz(
    symbol_letters="abc") == cb_symbol_abc
  #
  uhmx = "C 1 2 1 (%s)" % cb_symbol_xyz
  uhma = "C 1 2 1 (%s)" % cb_symbol_abc
  sx = sgtbx.space_group_info(symbol=uhmx)
  sa = sgtbx.space_group_info(symbol=uhma)
  assert sx.group() == sa.group()
  #
  # We trust that the cctbx is self-consistent.
  structure_unconv = random_structure.xray_structure(
    space_group_info=sx,
    elements=["C"],
    volume_per_atom=100,
    general_positions_only=True)
  assert str(structure_unconv.space_group_info()) == uhmx
  cb_op = structure_unconv.change_of_basis_op_to_reference_setting()
  structure_reference = structure_unconv.change_basis(cb_op=cb_op)
  assert str(structure_reference.space_group_info()) == "C 1 2 1"
  #
  # Verify the definitions in the paper based on the assumption
  # that the cctbx is self-consistent.
  site_reference = structure_reference.scatterers()[0].site
  site_unconv_direct = cb_mx_xyz * site_reference
  assert approx_equal(
    site_unconv_direct, structure_unconv.scatterers()[0].site)

def exercise_non_crystallographic_symmetry():
  sites_cart = flex.vec3_double(
    [(0.28730079491792732, 0.14711550696452974, 0.13031757579425293),
     (0.26144164573900441, 0.26385801128667269, 0.24113874888074088),
     (0.19728759424697784, 0.93346148983888833, 0.91783953828686837)])
  n = crystal.non_crystallographic_symmetry(sites_cart=sites_cart)
  assert approx_equal(n.unit_cell().parameters(),
    (1.6650571, 2.3613899, 2.36256589, 90, 90, 90))
  assert n.space_group_info().type().number() == 1
  n = crystal.non_crystallographic_symmetry(
    sites_cart=sites_cart, min_unit_cell_length=2)
  assert approx_equal(n.unit_cell().parameters(),
    (2, 2.3613899, 2.36256589, 90, 90, 90))
  sites_cart = flex.vec3_double(
    [(0.28730079491792732, 0.14711550696452974, 0.13031757579425293)])
  n = crystal.non_crystallographic_symmetry(sites_cart=sites_cart)
  assert approx_equal(n.unit_cell().parameters(), (1, 1, 1, 90, 90, 90))
  n = crystal.non_crystallographic_symmetry(
    sites_cart=sites_cart, default_buffer_layer=1.5)
  assert approx_equal(n.unit_cell().parameters(), (3, 3, 3, 90, 90, 90))

def exercise_special_position_settings():
  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  sp = crystal.special_position_settings(xs, 1, 2, False)
  assert sp.min_distance_sym_equiv() == 1
  assert sp.u_star_tolerance() == 2
  assert sp.assert_min_distance_sym_equiv() == False
  assert sp.site_symmetry((0,0,0)).multiplicity() == 1
  assert sp.site_symmetry(site=(0,0,0)).multiplicity() == 1
  assert sp.site_symmetry(site_cart=(0,0,0)).multiplicity() == 1
  assert str(sp.sym_equiv_sites((0,0,0)).special_op()) == "0,0,0"
  sites_cart = flex.vec3_double([(2,1,3), (0,0,0)])
  t = sp.site_symmetry_table(sites_cart=sites_cart)
  assert list(t.special_position_indices()) == [1]
  assert approx_equal(
    t.apply_symmetry_sites(unit_cell=xs.unit_cell(), sites_cart=sites_cart),
    sites_cart)
  #
  for min_distance_sym_equiv,special_op in [(1e-6, "0,0,0"), (0, "x,y,z")]:
    sp = crystal.special_position_settings(
      crystal_symmetry=xs,
      min_distance_sym_equiv=min_distance_sym_equiv)
    assert str(sp.sym_equiv_sites((0,0,0)).special_op()) == special_op
  #
  sites_cart = flex.vec3_double([(0,0,0)])
  sp = xs.special_position_settings()
  asu_mappings = sp.asu_mappings(buffer_thickness=3, sites_cart=sites_cart)
  assert list(asu_mappings.site_symmetry_table().special_position_indices()) \
      == [0]
  #
  pair_generator = sp.pair_generator(distance_cutoff=1, sites_cart=sites_cart)
  assert pair_generator.count_pairs() == 0
  sp0 = xs.special_position_settings(min_distance_sym_equiv=0)
  pair_generator = sp0.pair_generator(distance_cutoff=1, sites_cart=sites_cart)
  assert pair_generator.count_pairs() == 3
  #
  pair_asu_table = sp.pair_asu_table(distance_cutoff=1, sites_cart=sites_cart)
  assert pair_asu_table.table()[0].size() == 0
  pair_asu_table = sp0.pair_asu_table(distance_cutoff=1, sites_cart=sites_cart)
  assert pair_asu_table.table()[0][0].size() == 3

def exercise_site_symmetry(space_group_info):
  special_position_settings = crystal.special_position_settings(
    crystal_symmetry=space_group_info.any_compatible_crystal_symmetry(
      volume=1000))
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

def exercise_subtract_continuous_allowed_origin_shifts(
      space_group_info,
      use_niggli_cell,
      n_elements=3):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=300,
    min_distance=3.,
    general_positions_only=False)
  if (use_niggli_cell):
    structure = structure.niggli_cell()
  f_obs = abs(structure.structure_factors(
    d_min=3, algorithm="direct").f_calc())
  assert f_obs.indices().size() >= 10
  transl = matrix.col(flex.random_double_point_on_sphere()) * 2.345
  transl_no_cont = matrix.col(
    structure.subtract_continuous_allowed_origin_shifts(
      translation_cart=transl))
  transl_cont = transl - transl_no_cont
  structure_transl = structure.apply_shift(
    shift=structure.unit_cell().fractionalize(transl_cont),
    recompute_site_symmetries=True)
  f_transl = abs(f_obs.structure_factors_from_scatterers(
    xray_structure=structure_transl, algorithm="direct").f_calc())
  assert approx_equal(f_transl.data(), f_obs.data())

def run_call_back(flags, space_group_info):
  exercise_site_symmetry(space_group_info)
  for use_niggli_cell in [False, True]:
    exercise_subtract_continuous_allowed_origin_shifts(
      space_group_info=space_group_info,
      use_niggli_cell=use_niggli_cell)

def run():
  exercise_symmetry()
  exercise_correct_rhombohedral_setting_if_necessary()
  exercise_non_crystallographic_symmetry()
  exercise_special_position_settings()
  exercise_select_crystal_symmetry()
  verify_definitions_in_paper_zwart_2007()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
