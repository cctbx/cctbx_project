from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.model
from cctbx import geometry_restraints
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import iotbx
import sys
from cctbx.geometry_restraints.linking_class import linking_class
origin_ids = linking_class()


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
    dihedral {
      action = *add delete change
      atom_selection_1 = name C
      atom_selection_2 = name CA
      atom_selection_3 = name CB
      atom_selection_4 = name CG
      angle_ideal = 90
      sigma = 10
      periodicity = 1
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
  # Check the .geo file
  geo_fname = "pdb_interpretation_tst_edits_exercise_user_edits.geo"
  geometry.write_geo_file(file_name=geo_fname, site_labels=xrs.scatterers().extract_labels())
  with open(geo_fname, 'r') as f:
    user_suppl_count = 0
    for l in f.readlines():
      if l.find("| User supplied |")>-1:
        user_suppl_count += 1
    # Right now user-supplied planarity is missing from the .geo file.
    assert user_suppl_count == 5, "Expected 5 user-supplied restraints, got %i in the .geo" % user_suppl_count

  # initial
  assert geometry.pair_proxies().bond_proxies.simple.size() == 9
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.angle_proxies.size() == 9
  assert geometry.dihedral_proxies.size() == 7
  assert geometry.planarity_proxies.size() == 1
  assert geometry.parallelity_proxies.size() == 1

  ubonds_simpe, ubonds_asu, uangles, udihedrals, uplanarity, uparallelity = \
      geometry.get_user_supplied_restraints()
  assert ubonds_simpe.size() == 1
  assert ubonds_asu.size() == 0
  assert uangles.size() == 1
  assert udihedrals.size() == 1
  assert uplanarity.size() == 1
  assert uparallelity.size() == 1
  # make sure geometry stays the same
  assert geometry.pair_proxies().bond_proxies.simple.size() == 9
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.angle_proxies.size() == 9
  assert geometry.dihedral_proxies.size() == 7
  assert geometry.planarity_proxies.size() == 1
  assert geometry.parallelity_proxies.size() == 1
  # test functions one by one
  simple, asu = geometry.get_bond_proxies_without_user_supplied()
  assert simple.size() == 8
  assert asu.size() == 0
  angle = geometry.get_angle_proxies_without_user_supplied()
  assert angle.size() == 8
  dihed = geometry.get_dihedral_proxies_without_user_supplied()
  assert dihed.size() == 6
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

raw_records2 = """\
CRYST1   17.963   15.643   19.171  90.00  90.00  90.00 P 1
SCALE1      0.055670  0.000000  0.000000        0.00000
SCALE2      0.000000  0.063926  0.000000        0.00000
SCALE3      0.000000  0.000000  0.052162        0.00000
HETATM    1  N   ALA A   1      12.431   5.924  12.511  1.00 20.00      A    N
HETATM    2  CA  ALA A   1      11.018   6.145  12.230  1.00 20.00      A    C
HETATM    3  C   ALA A   1      10.790   7.554  11.693  1.00 20.00      A    C
HETATM    4  O   ALA A   1      11.665   8.411  11.803  1.00 20.00      A    O
HETATM    5  CB  ALA A   1      10.187   5.914  13.482  1.00 20.00      A    C
HETATM   13  N   ALA A   2       9.597   7.776  11.142  1.00 20.00      A    N
HETATM   14  CA  ALA A   2       9.208   9.039  10.514  1.00 20.00      A    C
HETATM   15  C   ALA A   2       8.148   8.584   9.515  1.00 20.00      A    C
HETATM   16  O   ALA A   2       7.467   7.584   9.741  1.00 20.00      A    O
HETATM   17  CB  ALA A   2      10.376   9.748   9.832  1.00 20.00      A    C
HETATM   23  N   ALA A   3       8.007   9.313   8.411  1.00 20.00      A    N
HETATM   24  CA  ALA A   3       7.032   8.963   7.385  1.00 20.00      A    C
HETATM   25  C   ALA A   3       7.369   9.651   6.067  1.00 20.00      A    C
HETATM   26  O   ALA A   3       8.098  10.643   6.037  1.00 20.00      A    O
HETATM   27  CB  ALA A   3       5.630   9.338   7.836  1.00 20.00      A    C
HETATM   28  OXT ALA A   3       6.921   9.231   5.000  1.00 20.00      A    O1-
""".splitlines()

def exercise_angle_edits_change(mon_lib_srv, ener_lib):
  edits = """\
refinement.geometry_restraints.edits {
  n_2_selection = chain A and resname ALA and resid 2 and name N
  ca_2_selection = chain A and resname ALA and resid 2 and name CA
  c_2_selection = chain A and resname ALA and resid 2 and name C
  angle {
    action = *change
    atom_selection_1 = $n_2_selection
    atom_selection_2 = $ca_2_selection
    atom_selection_3 = $c_2_selection
    angle_ideal = 100.00
    sigma = 5
  }
}"""
  gm_phil = iotbx.phil.parse(
      monomer_library.pdb_interpretation.grand_master_phil_str,
      process_includes=True)
  edits_phil = iotbx.phil.parse(edits)
  working_phil = gm_phil.fetch(edits_phil)
  params = working_phil.extract()
  # print params.geometry_restraints.edits.parallelity[0].atom_selection_1
  assert params.geometry_restraints.edits.angle[0].atom_selection_1.find("name N")
  processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records2,
      params = params.pdb_interpretation,
      log=None)
  grm = processed_pdb_file.geometry_restraints_manager(
      params_edits=params.geometry_restraints.edits,
      params_remove=params.geometry_restraints.remove)
  assert grm.angle_proxies.size() == 20
  user_defined = grm.angle_proxies.proxy_select(origin_id=origin_ids.get_origin_id('edits'))
  assert user_defined.size() == 1
  udp = user_defined[0]
  assert list(udp.i_seqs) == [5,6,7]
  assert approx_equal(udp.angle_ideal, 100, eps=1e-4)
  assert approx_equal(udp.weight, 0.04, eps=1e-4)

  from libtbx.test_utils import open_tmp_file
  from libtbx import easy_run
  pdb_file = open_tmp_file(suffix="aaa.pdb")
  pdb_file.write('\n'.join(raw_records2))
  pdb_file.close()
  edits_file = open_tmp_file(suffix="tau.edits")
  edits_file.write(edits)
  edits_file.close()
  cmd = "phenix.pdb_interpretation \"%s\" \"%s\" write_geo_files=True" %(
    pdb_file.name, edits_file.name)
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  geo_file = open(pdb_file.name+'.geo', "r")
  # geo_file = open(pdb_file.name.replace(".pdb", '_minimized.geo'), "r")
  geo_file_str = geo_file.read()
  assert '''| User supplied | restraints: 1
Sorted by residual:
angle pdb=" N   ALA A   2 " segid="A   "
      pdb=" CA  ALA A   2 " segid="A   "
      pdb=" C   ALA A   2 " segid="A   "
    ideal   model   delta    sigma   weight residual
   100.00''' in geo_file_str

def exercise_dihedral_edits_change():
  edits = """\
geometry_restraints.edits {
  dihedral {
    action = *change
    atom_selection_1 = resid 1 and name CA
    atom_selection_2 = resid 1 and name C
    atom_selection_3 = resid 2 and name N
    atom_selection_4 = resid 2 and name CA
    angle_ideal = 100.00
    alt_angle_ideals = 90,110
    sigma = 2
    periodicity = 2
  }
}"""
  def_params = mmtbx.model.manager.get_default_pdb_interpretation_scope()
  edits_phil = iotbx.phil.parse(edits)
  working_phil = def_params.fetch(edits_phil)
  params = working_phil.extract()
  inp = iotbx.pdb.input(lines=raw_records2, source_info=None)
  model = mmtbx.model.manager(model_input = inp)
  model.process(make_restraints=True)
  grm = model.get_restraints_manager().geometry
  assert grm.dihedral_proxies.size() == 9
  dih1 = grm.dihedral_proxies.proxy_select(
      n_seq=model.get_number_of_atoms(),
      iselection=flex.size_t([1,2,5,6]))
  assert dih1.size() == 1
  dp1 = dih1[0]
  assert approx_equal(dp1.angle_ideal, 180)
  assert dp1.alt_angle_ideals == None
  assert dp1.origin_id == 0
  assert dp1.periodicity == 0
  assert dp1.slack == 0
  assert not dp1.top_out
  assert approx_equal(dp1.weight, 0.04)

  # Now with modifications
  inp = iotbx.pdb.input(lines=raw_records2, source_info=None)
  model2 = mmtbx.model.manager(model_input = inp)
  model2.process(pdb_interpretation_params=params,
    make_restraints=True)
  grm2 = model2.get_restraints_manager().geometry
  assert grm2.dihedral_proxies.size() == 9
  dih2 = grm2.dihedral_proxies.proxy_select(
      n_seq=model2.get_number_of_atoms(),
      iselection=flex.size_t([1,2,5,6]))
  assert dih2.size() == 1
  dp2 = dih2[0]
  assert approx_equal(dp2.angle_ideal, 100)
  assert dp2.alt_angle_ideals == (90,110)
  assert dp2.origin_id == origin_ids.get_origin_id('edits')
  assert dp2.periodicity == 2
  assert dp2.slack == 0
  assert not dp2.top_out
  assert approx_equal(dp2.weight, 0.25)

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
    exercise_angle_edits_change(mon_lib_srv, ener_lib)
    exercise_dihedral_edits_change()

if (__name__ == "__main__"):
  exercise()
  print("OK")
