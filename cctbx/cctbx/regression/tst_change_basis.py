from cctbx.xray import minimization
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
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

def exercise(
      space_group_info,
      anomalous_flag,
      anisotropic_flag,
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
    anisotropic_flag=anisotropic_flag,
    random_occupancy=True)
  check_weight_without_occupancy(structure_z)
  f_abs_z = abs(structure_z.structure_factors(
    anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct").f_calc())
  if (0 or verbose):
    structure_z.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_z.special_position_indices().size()
  z2p_op = structure_z.space_group().z2p_op()
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
  f_abs_pz = abs(f_abs_z.structure_factors_from_scatterers(
    xray_structure=structure_pz,
    algorithm="direct").f_calc())
  c = flex.linear_correlation(f_abs_z.data(), f_abs_pz.data())
  assert c.is_well_defined()
  if (0 or verbose):
    print "correlation:", c.coefficient()
  assert c.coefficient() > 0.999
  f_abs_p_cb = f_abs_z.change_basis(z2p_op)
  o = flex.order(f_abs_p_cb.indices(), f_abs_z.indices())
  if (f_abs_z.space_group().n_ltr() == 1):
    assert o == 0
  else:
    assert o != 0
  f_abs_pz = f_abs_p_cb.change_basis(z2p_op.inverse())
  assert flex.order(f_abs_pz.indices(), f_abs_z.indices()) == 0
  f_abs_p_sf = abs(f_abs_p_cb.structure_factors_from_scatterers(
    xray_structure=structure_p,
    algorithm="direct").f_calc())
  c = flex.linear_correlation(f_abs_p_sf.data(), f_abs_p_cb.data())
  assert c.is_well_defined()
  if (0 or verbose):
    print "correlation:", c.coefficient()
  assert c.coefficient() > 0.999

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True)[:]: #SWITCH
    for anisotropic_flag in (False, True)[:]: #SWITCH
      exercise(
        space_group_info=space_group_info,
        anomalous_flag=anomalous_flag,
        anisotropic_flag=anisotropic_flag,
        verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
