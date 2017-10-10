from __future__ import division
from __future__ import print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
#import libtbx.load_env
from cctbx import geometry_restraints
import iotbx
import sys

raw_records1 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016447  0.009496  0.000000        0.00000
SCALE2      0.000000  0.018992  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010309        0.00000
ATOM   1050  N   LYS A 135      31.992  14.930  -7.233  1.00  9.47           N
ATOM   1051  CA  LYS A 135      31.388  16.216  -7.637  1.00 12.89           C
ATOM   1052  C   LYS A 135      30.807  16.840  -6.406  1.00  6.47           C
ATOM   1053  O   LYS A 135      29.583  16.869  -6.191  1.00 15.74           O
ATOM   1054  CB  LYS A 135      30.263  16.059  -8.655  1.00 13.51           C
ATOM   1055  CG  LYS A 135      30.742  15.277  -9.843  1.00 16.23           C
ATOM   1056  CD  LYS A 135      29.612  15.131 -10.835  1.00 28.55           C
ATOM   1057  CE  LYS A 135      30.173  14.812 -12.216  1.00 34.52           C
ATOM   1058  NZ  LYS A 135      29.396  13.756 -12.899  1.00 46.18           N
TER    1294      LYS A 162
END

""".splitlines()

def make_initial_grm(mon_lib_srv, ener_lib, records, params_edits):
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    # params         = params,
    raw_records    = records,
    force_symmetry = True,
    # log = sys.stdout,
    )

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = True,
    params_edits       = params_edits,
    plain_pairs_radius = 5.0,
    )
  xrs = processed_pdb_file.xray_structure()
  return geometry, xrs

def exercise_user_edits(mon_lib_srv, ener_lib):
  """
  exercise functions for managing geometry_restraints.edits scope for
  user-supplied restraints.
  """
  params_edits_str = """
    excessive_bond_distance_limit = 10
    bond {
      action = *add delete change
      atom_selection_1 = name CB
      atom_selection_2 = name CE
      symmetry_operation = None
      distance_ideal = 4
      sigma = 0.2
      slack = None
    }
    angle {
      action = *add delete change
      atom_selection_1 = name CB
      atom_selection_2 = name CD
      atom_selection_3 = name NZ
      angle_ideal = 120
      sigma = 3
    }
    planarity {
      action = *add delete change
      atom_selection = name CB or name CD or name CE
      sigma = 0.2
    }
    parallelity {
      action = *add delete change
      atom_selection_1 = name N or name CA or name O
      atom_selection_2 = name CB or name CD or name NZ
      sigma = 0.027
      target_angle_deg = 0
    }  """
  master_phil = iotbx.phil.parse(
      mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str)
  user_phil = iotbx.phil.parse(params_edits_str)
  working_phil = master_phil.fetch(sources=[user_phil])
  # working_phil.show()

  wp_extract = working_phil.extract()
  # print wp_extract.bond[0].atom_selection_1
  # STOP()
  # edits = None
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1, wp_extract)
  # initial
  assert geometry.pair_proxies().bond_proxies.simple.size() == 9
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.angle_proxies.size() == 9
  assert geometry.planarity_proxies.size() == 1
  assert geometry.parallelity_proxies.size() == 1

  ubonds_simpe, ubonds_asu, uangles, uplanarity, uparallelity = geometry.get_user_supplied_restraints()
  assert ubonds_simpe.size() == 1
  assert ubonds_asu.size() == 0
  assert uangles.size() == 1
  assert uplanarity.size() == 1
  assert uparallelity.size() == 1
  # make sure geometry stays the same
  assert geometry.pair_proxies().bond_proxies.simple.size() == 9
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.angle_proxies.size() == 9
  assert geometry.planarity_proxies.size() == 1
  assert geometry.parallelity_proxies.size() == 1
  # test functions one by one
  simple, asu = geometry.get_bond_proxies_without_user_supplied()
  assert simple.size() == 8
  assert asu.size() == 0
  angle = geometry.get_angle_proxies_without_user_supplied()
  assert angle.size() == 8
  plan = geometry.get_planarity_proxies_without_user_supplied()
  assert plan.size() == 0
  par = geometry.get_parallelity_proxies_without_user_supplied()
  assert par.size() == 0
  # make sure geometry stays the same
  assert geometry.pair_proxies().bond_proxies.simple.size() == 9
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.angle_proxies.size() == 9
  assert geometry.planarity_proxies.size() == 1
  assert geometry.parallelity_proxies.size() == 1


def exercise():
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None and ener_lib is not None:
    exercise_user_edits(mon_lib_srv, ener_lib)
if (__name__ == "__main__"):
  exercise()
  print("OK")