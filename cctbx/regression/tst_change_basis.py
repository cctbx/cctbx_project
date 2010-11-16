from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from cctbx.regression.tst_miller import generate_random_hl
import scitbx.math
from libtbx.test_utils import approx_equal
from stdlib import math
import sys

def check_weight_without_occupancy(structure):
  order_z = structure.space_group().order_z()
  for sc in structure.scatterers():
    assert approx_equal(
      sc.weight_without_occupancy(),
      float(sc.multiplicity()) / order_z)

def check_site_symmetry_table(structure_a, cb_op, structure_b):
  tab_a = structure_a.site_symmetry_table()
  tab_b = structure_b.site_symmetry_table()
  assert tab_a.indices().all_eq(tab_b.indices())
  assert tab_a.special_position_indices().all_eq(
         tab_b.special_position_indices())
  for oa,ob in zip(tab_a.table(),tab_b.table()):
    assert str(cb_op.apply(oa.special_op())) == str(ob.special_op())
    for ma,mb in zip(oa.matrices(),ob.matrices()):
      assert str(cb_op.apply(ma)) == str(mb)
    assert ob.multiplicity() > 0

def exercise(
      space_group_info,
      anomalous_flag,
      use_u_aniso,
      n_elements=3,
      d_min=3.,
      verbose=0):
  structure_z = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=200,
    random_f_prime_d_min=1.0,
    random_f_double_prime=anomalous_flag,
    random_u_iso=True,
    use_u_aniso=use_u_aniso,
    random_occupancy=True)
  check_weight_without_occupancy(structure_z)
  f_z = structure_z.structure_factors(
    anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct").f_calc()
  f_abs_z = abs(f_z)
  f_rad_z = f_z.phases()
  f_deg_z = f_z.phases(deg=True)
  hl_z = generate_random_hl(miller_set=f_z)
  hl_z_rad = hl_z.phase_integrals()
  if (0 or verbose):
    structure_z.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_z.special_position_indices().size()
  z2p_op = structure_z.space_group().z2p_op()
  z2p_op = sgtbx.change_of_basis_op(
      z2p_op.c()
    + sgtbx.tr_vec((2,-1,3), 12).new_denominator(z2p_op.c().t().den()))
  for change_hand in [False, True]:
    if (change_hand):
      z2p_op = z2p_op * sgtbx.change_of_basis_op("-x,-y,-z")
    structure_p = structure_z.change_basis(z2p_op)
    check_weight_without_occupancy(structure_p)
    check_site_symmetry_table(structure_z, z2p_op, structure_p)
    if (0 or verbose):
      structure_p.show_summary().show_scatterers()
      print "n_special_positions:", \
            structure_p.special_position_indices().size()
    assert tuple(structure_p.special_position_indices()) \
        == tuple(structure_z.special_position_indices())
    structure_pz = structure_p.change_basis(z2p_op.inverse())
    check_weight_without_occupancy(structure_pz)
    check_site_symmetry_table(structure_p, z2p_op.inverse(), structure_pz)
    assert structure_pz.unit_cell().is_similar_to(structure_z.unit_cell())
    assert structure_pz.space_group() == structure_z.space_group()
    f_pz = f_z.structure_factors_from_scatterers(
      xray_structure=structure_pz,
      algorithm="direct").f_calc()
    f_abs_pz = abs(f_pz)
    f_rad_pz = f_pz.phases()
    f_deg_pz = f_pz.phases(deg=True)
    c = flex.linear_correlation(f_abs_z.data(), f_abs_pz.data())
    assert c.is_well_defined()
    if (0 or verbose):
      print "correlation:", c.coefficient()
    assert c.coefficient() > 0.999
    f_p_cb = f_z.change_basis(z2p_op)
    f_abs_p_cb = f_abs_z.change_basis(z2p_op)
    f_rad_p_cb = f_rad_z.change_basis(z2p_op, deg=False)
    f_deg_p_cb = f_deg_z.change_basis(z2p_op, deg=True)
    hl_p_cb = hl_z.change_basis(z2p_op)
    hl_p_cb_rad = hl_p_cb.phase_integrals()
    assert approx_equal(
      hl_z_rad.change_basis(z2p_op).data(), hl_p_cb_rad.data())
    assert f_abs_p_cb.indices().all_eq(f_p_cb.indices())
    if (not change_hand):
      o = flex.order(f_abs_p_cb.indices(), f_abs_z.indices())
    else:
      o = flex.order(-f_abs_p_cb.indices(), f_abs_z.indices())
    if (f_abs_z.space_group().n_ltr() == 1):
      assert o == 0
    else:
      assert o != 0
    f_pz = f_p_cb.change_basis(z2p_op.inverse())
    f_abs_pz = f_abs_p_cb.change_basis(z2p_op.inverse())
    f_rad_pz = f_rad_p_cb.change_basis(z2p_op.inverse(), deg=False)
    f_deg_pz = f_deg_p_cb.change_basis(z2p_op.inverse(), deg=True)
    hl_pz = hl_p_cb.change_basis(z2p_op.inverse())
    hl_pz_rad = hl_pz.phase_integrals()
    assert approx_equal(hl_z_rad.data(), hl_pz_rad.data())
    for i,o in zip(hl_z.data(), hl_pz.data()):
      assert approx_equal(i, o)
    assert approx_equal(
      flex.max(flex.abs(f_pz.data() - f_z.data())), 0)
    assert flex.order(f_abs_pz.indices(), f_abs_z.indices()) == 0
    assert f_abs_pz.indices().all_eq(f_pz.indices())
    assert approx_equal(flex.max(scitbx.math.phase_error(
      phi1=f_rad_pz.data(), phi2=f_rad_z.data())), 0)
    assert approx_equal(flex.max(scitbx.math.phase_error(
      phi1=f_deg_pz.data(), phi2=f_deg_z.data(), deg=True)), 0)
    assert approx_equal(f_deg_pz.data(), f_rad_pz.data()*(180/math.pi))
    f_p_sf = f_p_cb.structure_factors_from_scatterers(
      xray_structure=structure_p,
      algorithm="direct").f_calc()
    det = abs(z2p_op.c().r().determinant())
    assert approx_equal(
      flex.max(flex.abs(f_p_sf.data()*complex(det) - f_p_cb.data())), 0)
    f_abs_p_sf = abs(f_p_sf)
    f_rad_p_sf = f_p_sf.phases()
    f_deg_p_sf = f_p_sf.phases(deg=True)
    assert approx_equal(flex.max(scitbx.math.phase_error(
      phi1=f_rad_p_sf.data(), phi2=f_rad_p_cb.data())), 0)
    assert approx_equal(flex.max(scitbx.math.phase_error(
      phi1=f_deg_p_sf.data(), phi2=f_deg_p_cb.data(), deg=True)), 0)
    c = flex.linear_correlation(f_abs_p_sf.data(), f_abs_p_cb.data())
    assert c.is_well_defined()
    if (0 or verbose):
      print "correlation:", c.coefficient()
    assert c.coefficient() > 0.999

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True)[:]: #SWITCH
    for use_u_aniso in (False, True)[:]: #SWITCH
      exercise(
        space_group_info=space_group_info,
        anomalous_flag=anomalous_flag,
        use_u_aniso=use_u_aniso,
        verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
