from cctbx.xray import minimization
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys

def exercise(space_group_info, anomalous_flag,
             n_elements=3, d_min=3., verbose=0):
  structure_z = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=200,
    random_f_prime_d_min=1.0,
    random_f_double_prime=anomalous_flag,
    random_u_iso=0001,
    random_occupancy=0001)
  f_abs_z_array = abs(structure_z.structure_factors_direct(
    anomalous_flag=anomalous_flag, d_min=d_min).f_calc_array())
  if (0 or verbose):
    structure_z.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_z.special_position_indices().size()
  z2p_op = structure_z.space_group().z2p_op()
  structure_p = structure_z.change_basis(z2p_op)
  if (0 or verbose):
    structure_p.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_p.special_position_indices().size()
  assert tuple(structure_p.special_position_indices()) \
      == tuple(structure_z.special_position_indices())
  structure_pz = structure_p.change_basis(z2p_op.inverse())
  assert structure_pz.unit_cell().is_similar_to(structure_z.unit_cell())
  assert structure_pz.space_group() == structure_z.space_group()
  f_abs_pz_array = abs(xray.structure_factors_direct(
    xray_structure=structure_pz,
    miller_set=f_abs_z_array).f_calc_array())
  r = flex.linear_regression(f_abs_z_array.data(), f_abs_pz_array.data())
  assert r.is_well_defined()
  if (0 or verbose):
    print "correlation:", r.correlation()
  assert r.correlation() > 0.999
  f_abs_p_array_cb = f_abs_z_array.change_basis(z2p_op)
  o = flex.order(f_abs_p_array_cb.indices(), f_abs_z_array.indices())
  if (f_abs_z_array.space_group().n_ltr() == 1):
    assert o == 0
  else:
    assert o != 0
  f_abs_pz_array = f_abs_p_array_cb.change_basis(z2p_op.inverse())
  assert flex.order(f_abs_pz_array.indices(), f_abs_z_array.indices()) == 0
  f_abs_p_array_sf = abs(xray.structure_factors_direct(
    xray_structure=structure_p,
    miller_set=f_abs_p_array_cb).f_calc_array())
  r = flex.linear_regression(f_abs_p_array_sf.data(), f_abs_p_array_cb.data())
  assert r.is_well_defined()
  if (0 or verbose):
    print "correlation:", r.correlation()
  assert r.correlation() > 0.999

def run_call_back(flags, space_group_info):
  for anomalous_flag in (00000, 0001)[:]: #SWITCH
    exercise(space_group_info, anomalous_flag, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
