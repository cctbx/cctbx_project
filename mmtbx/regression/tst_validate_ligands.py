from __future__ import absolute_import, division, print_function
import time, os, traceback
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
import libtbx.load_env
from libtbx.test_utils import approx_equal
from iotbx.cli_parser import run_program
from mmtbx.programs import validate_ligands as val_lig

# ------------------------------------------------------------------------------

def run():
  run_test1()
  run_test2()
  run_test3()
  #run_test4()
  run_test5()
  run_test6()

# ------------------------------------------------------------------------------

def run_test1():
  '''
  Several tests:
    - check if iselection for ligand PG5 (chain A resseq 201) is correct
    - count clashes involving ligand and calculate ligand clashscore
    - check geometry outliers
  '''
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/one_chain_ligand_water.pdb",
    test=os.path.isfile)
  args=[pdb_fname]
  print("mmtbx.development.validate_ligands %s" %(" ".join(args)))
  try:
    result = run_program(program_class=val_lig.Program,args=args,
     logger = null_out())
  except Exception as e:
    msg = traceback.format_exc()

  vl_manager = result.ligand_manager

  assert (len(vl_manager) == 1)
  lr = vl_manager[0]
  assert (lr.id_str == 'PG5 A 201')

  # test iselection
  assert(list(lr.ligand_isel) == [190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
    200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214,
    215, 216, 217, 218, 219])

  # Number of clashes and clashscore
  clashes_result = lr.get_overlaps()
  assert(clashes_result.n_clashes == 5)
  assert approx_equal(clashes_result.clashscore, 27.6, eps=0.5)
  #
  # Geometry deviations
  rmsd_result = lr.get_rmsds()
  assert approx_equal(rmsd_result.bond_rmsd, 0.093, eps=0.005)
  assert approx_equal(rmsd_result.bond_rmsz, 4.654, eps=0.005)
  assert(rmsd_result.bond_n_outliers == 4)
  assert approx_equal(rmsd_result.angle_rmsd, 4.48, eps=0.05)
  assert approx_equal(rmsd_result.angle_rmsz, 1.49, eps=0.05)
  assert(rmsd_result.angle_n_outliers == 0)
  assert approx_equal(rmsd_result.dihedral_rmsd, 32.6, eps=0.5)

  os.remove('one_chain_ligand_water_newH.cif')
  os.remove('one_chain_ligand_water_newH.txt')

# ------------------------------------------------------------------------------

def run_test2():
  '''
  Test
  - occupancy determination for ligands
  - adp determination for ligands and neighbors
  '''
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/two_chains_ligand_water.pdb",
    test=os.path.isfile)
  args=[pdb_fname]
  print("mmtbx.development.validate_ligands %s" %(" ".join(args)))
  try:
    result = run_program(program_class=val_lig.Program,args=args,
     logger = null_out())
  except Exception as e:
    msg = traceback.format_exc()

  vl_manager = result.ligand_manager
  tst_occupancies(vl_manager = vl_manager)
  tst_adps(vl_manager = vl_manager)
  tst_rmsds(vl_manager = vl_manager)

  os.remove('two_chains_ligand_water_newH.cif')
  os.remove('two_chains_ligand_water_newH.txt')

# ------------------------------------------------------------------------------

def tst_rmsds(vl_manager):
  for lr in vl_manager:
    id_str = lr.id_str
    if (id_str.strip() == 'SO4 A   4'):
      rmsd_result = lr.get_rmsds()
      assert approx_equal(rmsd_result.bond_rmsd, 0.055, eps=0.005)
      assert approx_equal(rmsd_result.bond_rmsz, 2.761, eps=0.005)
      assert(rmsd_result.bond_n_outliers == 1)
      #
      assert approx_equal(rmsd_result.angle_rmsd, 4.08, eps=0.05)
      assert approx_equal(rmsd_result.angle_rmsz, 1.36, eps=0.05)
      assert(rmsd_result.angle_n_outliers == 0)

# ------------------------------------------------------------------------------

def tst_occupancies(vl_manager):
  '''
  Test occupancies
  '''
  assert (len(vl_manager) == 5)
  for lr in vl_manager:
    occs = lr.get_occupancies()
    id_str = lr.id_str
    if (id_str.strip() == 'ABEN A   2'):
      assert approx_equal(occs.occ_mean, 0.56, eps=0.01)
    if (id_str.strip() == 'BBEN A   2'):
      assert approx_equal(occs.occ_mean, 0.44, eps=0.01)
    if (id_str.strip() == 'SO4 A   3'):
      assert approx_equal(occs.occ_mean, 0.65, eps=0.01)
    if (id_str.strip() == 'SO4 A   4'):
      assert approx_equal(occs.occ_mean, 0.48, eps=0.01)
    if (id_str.strip() == 'GOL A   5'):
      assert approx_equal(occs.occ_mean, 0.67, eps=0.01)

# ------------------------------------------------------------------------------

def tst_adps(vl_manager):
  '''
  Test ADPs of ligands and surrounding atoms
  '''
  for lr in vl_manager:
    occs = lr.get_occupancies()
    id_str = lr.id_str
    adps = lr.get_adps()
    if (id_str == 'ABEN B   2'):
      assert(adps.n_iso == 0)
      assert(adps.n_aniso == 9)
      assert(adps.b_min_within is None)
      assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
        [4.7, 8.1, 6.0], eps=0.1)
    if (id_str.strip() == 'BBEN B   2'):
      assert(adps.n_iso == 0)
      assert(adps.n_aniso == 9)
      assert(adps.b_min_within is None)
      assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
        [5.1, 8.2, 6.4], eps=0.1)
    if (id_str.strip() == 'SO4 B   3'):
      assert(adps.n_iso == 5)
      assert(adps.n_aniso == 0)
      assert(adps.b_min_within is None)
      assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
        [7.4,13.1,10.2], eps=0.1)
    if (id_str.strip() == 'SO4 B   4'):
      assert(adps.n_iso == 0)
      assert(adps.n_aniso == 5)
      assert(adps.b_min_within is None)
      assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
        [10.3,14.6,12.3], eps=0.1)
    if (id_str.strip() == 'GOL B   5'):
      assert(adps.n_iso == 6)
      assert(adps.n_aniso == 0)
      assert(adps.b_min_within is None)
      assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
        [58.7,114.9,96.9], eps=0.1)

# ------------------------------------------------------------------------------

def run_test3():
  '''
  Test
  - CC calculation for three ligands
  - ADP calculations for ligands
  - occupancy calculations for ligands
  '''
  mtz_fname = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1avd.mtz",
    test=os.path.isfile)
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1avd.ent.gz",
    test=os.path.isfile)
  args=[pdb_fname, mtz_fname]
  #
  print("mmtbx.development.validate_ligands %s" %(" ".join(args)))
  try:
    result = run_program(program_class=val_lig.Program,args=args,
     logger = null_out())
  except Exception as e:
    msg = traceback.format_exc()

  vl_manager = result.ligand_manager

  for lr in vl_manager:
    occs = lr.get_occupancies()
    id_str = lr.id_str
    adps = lr.get_adps()
    ccs = lr.get_ccs()
    clashes_result = lr.get_overlaps()
    #
    if (id_str.strip() == 'NAG A 600'):
      assert approx_equal(ccs.cc_2fofc, 0.87, eps=0.03)
      #
      assert approx_equal(occs.occ_min, 0, eps=0.01)
      assert approx_equal(occs.occ_max, 1, eps=0.01)
      assert approx_equal(occs.occ_mean, 0.86, eps=0.01)
      #
      assert approx_equal(adps.b_min, 27.99, eps=0.01)
      assert approx_equal(adps.b_max, 90.00, eps=0.01)
      assert approx_equal(adps.b_mean, 70.71, eps=0.05)
      #
      assert approx_equal(adps.b_min_within, 5.23, eps=0.01)
      assert approx_equal(adps.b_max_within, 79.14, eps=0.01)
      assert approx_equal(adps.b_mean_within, 35.37, eps=0.02)
      #
      assert(clashes_result.n_clashes == 1)
      assert approx_equal(clashes_result.clashscore, 9.4, eps=0.5)
      #
    if (id_str.strip() == 'BTN A 400'):
      assert approx_equal(ccs.cc_2fofc, 0.94, eps=0.03)
      #
      assert approx_equal(occs.occ_mean, 1, eps=0.01)
      #
      assert approx_equal(adps.b_min, 4.00, eps=0.01)
      assert approx_equal(adps.b_max, 90.00, eps=0.01)
      assert approx_equal(adps.b_mean, 31.19, eps=0.05)
      #
      assert approx_equal(adps.b_min_within, 4.00, eps=0.01)
      assert approx_equal(adps.b_max_within, 54.65, eps=0.01)
      assert approx_equal(adps.b_mean_within, 23.23, eps=0.02)
      #
      assert(clashes_result.n_clashes == 4)
      assert approx_equal(clashes_result.clashscore, 13.0, eps=0.5)

      rmsd_result = lr.get_rmsds()
      assert approx_equal(rmsd_result.bond_rmsd, 0.033, eps=0.005)
      assert approx_equal(rmsd_result.bond_rmsz, 1.639, eps=0.005)
      assert(rmsd_result.bond_n_outliers == 1)
      #
      assert approx_equal(rmsd_result.angle_rmsd, 4.00, eps=0.05)
      assert approx_equal(rmsd_result.angle_rmsz, 1.33, eps=0.05)
      assert(rmsd_result.angle_n_outliers == 1)
      assert approx_equal(rmsd_result.dihedral_rmsd, 17.4, eps=0.5)
      assert approx_equal(rmsd_result.planarity_rmsd, 0.02, eps=0.05)
      #
    if (id_str.strip() == 'BTN B 401'):
      assert approx_equal(ccs.cc_2fofc, 0.95, eps=0.03)
      #
      assert approx_equal(occs.occ_mean, 1, eps=0.01)
      #
      assert approx_equal(adps.b_min, 4.00, eps=0.01)
      assert approx_equal(adps.b_max, 46.67, eps=0.01)
      assert approx_equal(adps.b_mean, 23.04, eps=0.05)
      #
      assert approx_equal(adps.b_min_within, 4.00, eps=0.01)
      assert approx_equal(adps.b_max_within, 75.42, eps=0.01)
      assert approx_equal(adps.b_mean_within, 28.16, eps=0.02)
      #
      assert(clashes_result.n_clashes == 6)
      assert approx_equal(clashes_result.clashscore, 18.5, eps=0.5)

      #print(round(clashes_result.clashscore,1))
      #print(adps.b_min_within)
      #print(adps.b_max_within)
      #print(adps.b_mean_within)
      #print(occs.occ_min)
      #print(occs.occ_max)
      #print(occs.occ_mean)
      #print(adps.b_min)
      #print(adps.b_max)
      #print(adps.b_mean)

  os.remove('pdb1avd_newH.txt')
  os.remove('pdb1avd_newH.cif')

#def tst_get_overlaps(vl_manager):
#  '''
#  Test nonbonded overlaps
#  '''
#  for id_tuple, ligand_dict in vl_manager.items():
#    for altloc, lr in ligand_dict.items():
#      clashes_result = lr.get_overlaps()
#      assert(clashes_result.n_clashes == 5)
#      assert approx_equal(clashes_result.clashscore, 31.6, eps=1.0)
# anaconda
#(['pdb=" HE3 MET A 107 "', 'pdb=" H81 PG5 A 201 "'], 17, 54, 2.0370952358689647, 2.44, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" C8  PG5 A 201 "'], 15, 34, 2.946989989803154, 3.4, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" H83 PG5 A 201 "'], 15, 56, 2.4839921497460486, 2.92, '', None)
#
# MAC
#(['pdb=" CE  MET A 107 "', 'pdb=" C8  PG5 A 201 "'], 16, 35, 2.946989989803154, 3.4, '', None),
#(['pdb=" HE3 MET A 107 "', 'pdb=" H83 PG5 A 201 "'], 18, 57, 2.026073542594147, 2.44, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" H81 PG5 A 201 "'], 16, 55, 2.4973179613337146, 2.92, '', None)
#
# Readyset gives different names to H atoms.

# ------------------------------------------------------------------------------

def run_test4():
  '''
  Test if ligands with 5-letter ID are processed properly
  from PDB model 7has
  '''
  model_fn = "tst_4_fragment.cif"
  with open(model_fn, "w") as f:
    f.write(cif_str_tst_4)
  args = [model_fn]
  print("mmtbx.development.validate_ligands %s run_reduce2=False" % model_fn)
  try:
    result = run_program(program_class=val_lig.Program,args=args,
     logger = null_out())
  except Exception as e:
    msg = traceback.format_exc()

  vl_manager = result.ligand_manager

  for lr in vl_manager:
    assert (lr.id_str == 'A1AYY A 301')
    occs = lr.get_occupancies()
    adps = lr.get_adps()
    #
    assert approx_equal(occs.occ_mean, 0.91, eps=0.01)
    #
    assert approx_equal(adps.b_min, 38.36, eps=0.01)
    assert approx_equal(adps.b_max, 66.40, eps=0.01)
    assert approx_equal(adps.b_mean, 47.82, eps=0.05)
    #
    assert approx_equal(adps.b_min_within, 30.55, eps=0.01)
    assert approx_equal(adps.b_max_within, 63.10, eps=0.01)
    assert approx_equal(adps.b_mean_within, 42.45, eps=0.02)
  #
  os.remove(model_fn)
  os.remove('tst_4_fragment_newH.cif')
  os.remove('tst_4_fragment_newH.txt')

# ------------------------------------------------------------------------------

def run_test5():
  '''
  Test if ligand with two conformers with different names is processed correctly
  from PDB model 4d3w
  '''
  model_fn = "tst_5.pdb"
  with open(model_fn, "w") as f:
    f.write(pdb_str_tst_5)
  args = [model_fn]
  print("mmtbx.development.validate_ligands %s" % model_fn)
  try:
    result = run_program(program_class=val_lig.Program,args=args,
     logger = null_out())
  except Exception as e:
    msg = traceback.format_exc()

  vl_manager = result.ligand_manager

  for lr in vl_manager:
    assert (lr.id_str in ['ANGA B   1', 'BA2G B   1', 'GAL B   2'])
  #
  os.remove(model_fn)
  os.remove('tst_5_newH.cif')
  os.remove('tst_5_newH.txt')

# ------------------------------------------------------------------------------

def run_test6():
  pdb_str = '''
CRYST1   14.103   13.596   12.544  90.00  90.00  90.00 P 1
SCALE1      0.070907  0.000000  0.000000        0.00000
SCALE2      0.000000  0.073551  0.000000        0.00000
SCALE3      0.000000  0.000000  0.079719        0.00000
ATOM      1  N   SER A   9       8.404   8.079   5.179  1.00 16.08           N
ATOM      2  CA  SER A   9       7.167   7.545   5.742  1.00 17.80           C
ATOM      3  C   SER A   9       6.169   8.596   6.220  1.00 19.90           C
ATOM      4  O   SER A   9       5.000   8.269   6.456  1.00 19.71           O
ATOM      5  CB  SER A   9       7.456   6.619   6.919  1.00 16.95           C
ATOM      6  OG  SER A   9       8.035   5.415   6.447  1.00 17.71           O
ATOM      7  H   SER A   9       9.103   7.949   5.663  1.00 16.08           H
ATOM      8  HA  SER A   9       6.761   7.070   5.000  1.00 17.80           H
ATOM      9  HB2 SER A   9       8.057   7.054   7.544  1.00 16.95           H
ATOM     10  HB3 SER A   9       6.636   6.427   7.400  1.00 16.95           H
ATOM     11  HG  SER A   9       8.382   5.000   7.090  1.00 17.71           H
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  stats = model.geometry_statistics(use_hydrogens=True)
  dihedral = stats.dihedral()
  dihedralz = stats.dihedral(return_rmsZ=True)
  assert dihedral.n == 3
  assert len(dihedral.outliers) == 0
  assert approx_equal(dihedral.mean, 10.91, eps=0.05)
  assert approx_equal(dihedralz.mean, 0.52, eps=0.05)

# ------------------------------------------------------------------------------

cif_str_tst_4 = '''
data_default
_cell.length_a                    67.330
_cell.length_b                    91.060
_cell.length_c                    98.090
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      601396.637
_space_group.crystal_system       orthorhombic
_space_group.IT_number            23
_space_group.name_H-M_alt         'I 2 2 2'
_space_group.name_Hall            ' I 2 2'
_symmetry.space_group_name_H-M    'I 2 2 2'
_symmetry.space_group_name_Hall   ' I 2 2'
_symmetry.Int_Tables_number       23
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z
   2 x,-y,-z
   3 -x,y,-z
   4 -x,-y,z
   5 x+1/2,y+1/2,z+1/2
   6 x+1/2,-y+1/2,-z+1/2
   7 -x+1/2,y+1/2,-z+1/2
   8 -x+1/2,-y+1/2,z+1/2

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   ATOM 1452 N . LEU A 107 ? 9.35882 -9.46004 -29.05718 1.000 40.04006 N ? A ? 92 N 1
   ATOM 1453 CA . LEU A 107 ? 9.05649 -9.42408 -27.62920 1.000 44.50060 C ? A ? 92 CA 1
   ATOM 1454 C . LEU A 107 ? 10.31164 -9.61463 -26.78677 1.000 46.34740 C ? A ? 92 C 1
   ATOM 1455 O . LEU A 107 ? 10.50736 -8.91194 -25.78853 1.000 46.42741 O ? A ? 92 O 1
   ATOM 1456 CB . LEU A 107 ? 8.01963 -10.49032 -27.29190 1.000 38.73807 C ? A ? 92 CB 1
   ATOM 1457 CG . LEU A 107 ? 6.62315 -10.23446 -27.85167 1.000 43.87342 C ? A ? 92 CG 1
   ATOM 1458 CD1 . LEU A 107 ? 5.73330 -11.41418 -27.53992 1.000 45.90445 C ? A ? 92 CD1 1
   ATOM 1459 CD2 . LEU A 107 ? 6.03694 -8.94954 -27.28925 1.000 39.54185 C ? A ? 92 CD2 1
   ATOM 1460 H . LEU A 107 ? 9.03211 -10.13489 -29.47842 1.000 40.04006 H ? A ? 92 H 1
   ATOM 1461 HA . LEU A 107 ? 8.68818 -8.55585 -27.40243 1.000 44.50060 H ? A ? 92 HA 1
   ATOM 1462 HB2 . LEU A 107 ? 7.94074 -10.54597 -26.32672 1.000 38.73807 H ? A ? 92 HB2 1
   ATOM 1463 HB3 . LEU A 107 ? 8.32562 -11.33861 -27.64921 1.000 38.73807 H ? A ? 92 HB3 1
   ATOM 1464 HG . LEU A 107 ? 6.67444 -10.12742 -28.81438 1.000 43.87342 H ? A ? 92 HG 1
   ATOM 1465 HD11 . LEU A 107 ? 5.68474 -11.52758 -26.57780 1.000 45.90445 H ? A ? 92 HD11 1
   ATOM 1466 HD12 . LEU A 107 ? 6.10907 -12.20979 -27.94821 1.000 45.90445 H ? A ? 92 HD12 1
   ATOM 1467 HD13 . LEU A 107 ? 4.84786 -11.24525 -27.89818 1.000 45.90445 H ? A ? 92 HD13 1
   ATOM 1468 HD21 . LEU A 107 ? 6.61336 -8.20750 -27.53013 1.000 39.54185 H ? A ? 92 HD21 1
   ATOM 1469 HD22 . LEU A 107 ? 5.98014 -9.02441 -26.32381 1.000 39.54185 H ? A ? 92 HD22 1
   ATOM 1470 HD23 . LEU A 107 ? 5.15199 -8.81593 -27.66329 1.000 39.54185 H ? A ? 92 HD23 1
   ATOM 1511 N . ALA A 111 ? 11.99545 -7.07173 -24.18626 1.000 41.28916 N ? A ? 96 N 1
   ATOM 1512 CA . ALA A 111 ? 12.14347 -7.55654 -22.81984 1.000 50.33630 C ? A ? 96 CA 1
   ATOM 1513 C . ALA A 111 ? 13.58717 -7.55437 -22.33204 1.000 56.57231 C ? A ? 96 C 1
   ATOM 1514 O . ALA A 111 ? 13.80955 -7.73541 -21.12990 1.000 50.29603 O ? A ? 96 O 1
   ATOM 1515 CB . ALA A 111 ? 11.56279 -8.96693 -22.70298 1.000 40.65961 C ? A ? 96 CB 1
   ATOM 1516 H . ALA A 111 ? 11.80242 -7.68526 -24.75715 1.000 41.28916 H ? A ? 96 H 1
   ATOM 1517 HA . ALA A 111 ? 11.65287 -6.95070 -22.24264 1.000 50.33630 H ? A ? 96 HA 1
   ATOM 1518 HB1 . ALA A 111 ? 12.03907 -9.55576 -23.30906 1.000 40.65961 H ? A ? 96 HB1 1
   ATOM 1519 HB2 . ALA A 111 ? 10.62194 -8.94153 -22.93762 1.000 40.65961 H ? A ? 96 HB2 1
   ATOM 1520 HB3 . ALA A 111 ? 11.66773 -9.27700 -21.78988 1.000 40.65961 H ? A ? 96 HB3 1
   ATOM 1904 N . TYR A 139 ? 6.63097 -13.56248 -17.86012 1.000 30.55156 N ? A ? 124 N 1
   ATOM 1905 CA . TYR A 139 ? 7.45336 -14.64062 -17.32367 1.000 31.38192 C ? A ? 124 CA 1
   ATOM 1906 C . TYR A 139 ? 7.10974 -14.98654 -15.88117 1.000 32.69682 C ? A ? 124 C 1
   ATOM 1907 O . TYR A 139 ? 7.55705 -16.02976 -15.39099 1.000 36.48358 O ? A ? 124 O 1
   ATOM 1908 CB . TYR A 139 ? 8.93112 -14.27595 -17.45657 1.000 32.07174 C ? A ? 124 CB 1
   ATOM 1909 CG . TYR A 139 ? 9.33068 -14.07862 -18.89825 1.000 38.57699 C ? A ? 124 CG 1
   ATOM 1910 CD1 . TYR A 139 ? 9.44155 -15.16247 -19.76172 1.000 40.13981 C ? A ? 124 CD1 1
   ATOM 1911 CD2 . TYR A 139 ? 9.56183 -12.80879 -19.40739 1.000 41.42076 C ? A ? 124 CD2 1
   ATOM 1912 CE1 . TYR A 139 ? 9.79161 -14.98764 -21.09024 1.000 45.33438 C ? A ? 124 CE1 1
   ATOM 1913 CE2 . TYR A 139 ? 9.91275 -12.62417 -20.73234 1.000 42.40061 C ? A ? 124 CE2 1
   ATOM 1914 CZ . TYR A 139 ? 10.02524 -13.71558 -21.56910 1.000 42.57300 C ? A ? 124 CZ 1
   ATOM 1915 OH . TYR A 139 ? 10.37367 -13.53282 -22.88819 1.000 63.10415 O ? A ? 124 OH 1
   ATOM 1916 H . TYR A 139 ? 7.04775 -12.82240 -17.99492 1.000 30.55156 H ? A ? 124 H 1
   ATOM 1917 HA . TYR A 139 ? 7.28145 -15.44460 -17.83842 1.000 31.38192 H ? A ? 124 HA 1
   ATOM 1918 HB2 . TYR A 139 ? 9.10075 -13.44967 -16.97762 1.000 32.07174 H ? A ? 124 HB2 1
   ATOM 1919 HB3 . TYR A 139 ? 9.47142 -14.99119 -17.08589 1.000 32.07174 H ? A ? 124 HB3 1
   ATOM 1920 HD1 . TYR A 139 ? 9.27775 -16.02013 -19.44158 1.000 40.13981 H ? A ? 124 HD1 1
   ATOM 1921 HD2 . TYR A 139 ? 9.47946 -12.06989 -18.84869 1.000 41.42076 H ? A ? 124 HD2 1
   ATOM 1922 HE1 . TYR A 139 ? 9.86854 -15.72241 -21.65512 1.000 45.33438 H ? A ? 124 HE1 1
   ATOM 1923 HE2 . TYR A 139 ? 10.07228 -11.76780 -21.05805 1.000 42.40061 H ? A ? 124 HE2 1
   ATOM 1924 HH . TYR A 139 ? 10.08188 -12.79454 -23.16278 1.000 63.10415 H ? A ? 124 HH 1
   HETATM 3297 C1 . A1AYY A 301 ? 8.22925 -10.43941 -22.78041 0.908 57.20233 C ? B ? . C1 1
   HETATM 3298 C10 . A1AYY A 301 ? 6.78537 -14.46083 -25.10107 0.908 38.35600 C ? B ? . C10 1
   HETATM 3299 C2 . A1AYY A 301 ? 6.84230 -10.87462 -23.16045 0.908 50.12524 C ? B ? . C2 1
   HETATM 3300 C3 . A1AYY A 301 ? 5.71051 -10.11856 -22.88689 0.908 51.79520 C ? B ? . C3 1
   HETATM 3301 C4 . A1AYY A 301 ? 4.47648 -10.59742 -23.27872 0.908 47.07763 C ? B ? . C4 1
   HETATM 3302 C5 . A1AYY A 301 ? 4.38335 -11.81161 -23.93830 0.908 50.92592 C ? B ? . C5 1
   HETATM 3303 C6 . A1AYY A 301 ? 5.57038 -12.49750 -24.17023 0.908 45.47959 C ? B ? . C6 1
   HETATM 3304 C7 . A1AYY A 301 ? 4.38671 -14.37271 -25.31214 0.908 45.16716 C ? B ? . C7 1
   HETATM 3305 C8 . A1AYY A 301 ? 4.92432 -15.75041 -25.62964 0.908 39.66964 C ? B ? . C8 1
   HETATM 3306 C9 . A1AYY A 301 ? 6.29365 -15.42865 -26.15268 0.908 40.01913 C ? B ? . C9 1
   HETATM 3307 N1 . A1AYY A 301 ? 5.57761 -13.69344 -24.81157 0.908 40.78911 N ? B ? . N1 1
   HETATM 3308 N2 . A1AYY A 301 ? 6.77696 -12.05090 -23.79944 0.908 48.68856 N ? B ? . N2 1
   HETATM 3309 O1 . A1AYY A 301 ? 9.18903 -10.97661 -23.66412 0.908 66.39976 O ? B ? . O1 1
'''

pdb_str_tst_5 = '''
CRYST1   67.843   26.790   30.439  90.00  90.00  90.00 P 1
SCALE1      0.014740  0.000000  0.000000        0.00000
SCALE2      0.000000  0.037327  0.000000        0.00000
SCALE3      0.000000  0.000000  0.032853        0.00000
ATOM      1  N   TYR A  76      58.963   6.159  18.630  1.00 13.67           N
ATOM      2  C   TYR A  76      57.652   6.230  16.561  1.00 13.45           C
ATOM      3  O   TYR A  76      58.062   5.249  15.941  1.00 14.53           O
ATOM      4  CA ATYR A  76      58.559   7.008  17.518  0.50 14.05           C
ATOM      5  CB ATYR A  76      59.750   7.628  16.752  0.50 14.89           C
ATOM      6  CG ATYR A  76      59.283   8.584  15.663  0.50 16.30           C
ATOM      7  CD1ATYR A  76      58.642   9.783  15.996  0.50 17.47           C
ATOM      8  CD2ATYR A  76      59.437   8.280  14.316  0.50 17.43           C
ATOM      9  CE1ATYR A  76      58.192  10.650  15.022  0.50 18.29           C
ATOM     10  CE2ATYR A  76      58.987   9.137  13.335  0.50 17.65           C
ATOM     11  CZ ATYR A  76      58.364  10.319  13.691  0.50 19.22           C
ATOM     12  OH ATYR A  76      57.897  11.200  12.733  0.50 21.32           O
ATOM     13  CA BTYR A  76      58.566   6.992  17.509  0.50 13.57           C
ATOM     14  CB BTYR A  76      59.788   7.499  16.746  0.50 13.80           C
ATOM     15  CG BTYR A  76      60.605   8.433  17.575  0.50 14.39           C
ATOM     16  CD1BTYR A  76      61.568   7.949  18.448  0.50 14.45           C
ATOM     17  CD2BTYR A  76      60.411   9.803  17.508  0.50 15.07           C
ATOM     18  CE1BTYR A  76      62.316   8.804  19.222  0.50 15.57           C
ATOM     19  CE2BTYR A  76      61.170  10.664  18.280  0.50 15.82           C
ATOM     20  CZ BTYR A  76      62.115  10.158  19.127  0.50 15.67           C
ATOM     21  OH BTYR A  76      62.843  10.998  19.920  0.50 17.94           O
ATOM     22  N   PRO A  77      56.379   6.620  16.465  1.00 12.19           N
ATOM     23  CA  PRO A  77      55.680   7.664  17.202  1.00 11.58           C
ATOM     24  C   PRO A  77      55.507   7.326  18.675  1.00 10.59           C
ATOM     25  O   PRO A  77      55.674   6.178  19.090  1.00 10.13           O
ATOM     26  CB  PRO A  77      54.306   7.694  16.515  1.00 12.59           C
ATOM     27  CG  PRO A  77      54.542   7.081  15.184  1.00 13.41           C
ATOM     28  CD  PRO A  77      55.502   5.992  15.467  1.00 13.11           C
ATOM     29  N   TRP A  79      52.719   7.480  20.378  1.00  9.15           N
ATOM     30  CA  TRP A  79      51.452   6.851  20.599  1.00  9.12           C
ATOM     31  C   TRP A  79      51.516   5.378  20.145  1.00  9.16           C
ATOM     32  O   TRP A  79      52.373   5.000  19.354  1.00  9.65           O
ATOM     33  CB  TRP A  79      50.308   7.586  19.873  1.00  9.06           C
ATOM     34  CG  TRP A  79      50.559   7.920  18.458  1.00  8.62           C
ATOM     35  CD1 TRP A  79      51.101   9.076  17.989  1.00  8.79           C
ATOM     36  CD2 TRP A  79      50.326   7.108  17.311  1.00  8.00           C
ATOM     37  NE1 TRP A  79      51.203   9.049  16.634  1.00  8.50           N
ATOM     38  CE2 TRP A  79      50.763   7.825  16.188  1.00  8.27           C
ATOM     39  CE3 TRP A  79      49.775   5.840  17.115  1.00  7.63           C
ATOM     40  CZ2 TRP A  79      50.649   7.307  14.900  1.00  8.35           C
ATOM     41  CZ3 TRP A  79      49.647   5.359  15.833  1.00  8.08           C
ATOM     42  CH2 TRP A  79      50.094   6.075  14.759  1.00  8.28           C
ATOM     43  N   GLY A 118      52.212  17.874  22.406  1.00 10.97           N
ATOM     44  CA  GLY A 118      52.287  16.536  21.821  1.00  9.61           C
ATOM     45  C   GLY A 118      52.843  15.513  22.767  1.00 10.10           C
ATOM     46  O   GLY A 118      53.685  15.834  23.623  1.00  9.35           O
ATOM     47  N   CYS A 119      52.411  14.273  22.618  1.00  9.67           N
ATOM     48  CA  CYS A 119      53.033  13.132  23.277  1.00 10.60           C
ATOM     49  C   CYS A 119      53.102  13.263  24.794  1.00 10.47           C
ATOM     50  O   CYS A 119      53.997  12.760  25.439  1.00 10.67           O
ATOM     51  CB  CYS A 119      54.437  12.899  22.674  1.00 10.72           C
ATOM     52  SG  CYS A 119      54.337  12.396  20.924  1.00 12.24           S
ATOM     53  N   THR A 154       9.203   8.604  11.635  1.00  8.67           N
ATOM     54  CA  THR A 154       9.431   7.472  10.794  1.00  8.42           C
ATOM     55  C   THR A 154      10.007   7.956   9.480  1.00  8.45           C
ATOM     56  O   THR A 154       9.446   8.840   8.822  1.00  9.72           O
ATOM     57  CB  THR A 154       8.116   6.726  10.521  1.00  7.82           C
ATOM     58  OG1 THR A 154       7.641   6.150  11.740  1.00  7.74           O
ATOM     59  CG2 THR A 154       8.300   5.606   9.484  1.00  8.33           C
ATOM     60  N   GLY A 155      11.111   7.363   9.053  1.00  7.99           N
ATOM     61  CA  GLY A 155      11.710   7.739   7.798  1.00  8.58           C
ATOM     62  C   GLY A 155      13.216   7.596   7.781  1.00  8.60           C
ATOM     63  O   GLY A 155      13.809   7.064   8.710  1.00  8.47           O
ATOM     64  N   TRP A 198      40.718   9.388  12.154  1.00  7.65           N
ATOM     65  CA  TRP A 198      41.594  10.496  11.835  1.00  7.72           C
ATOM     66  C   TRP A 198      42.823   9.995  11.061  1.00  8.32           C
ATOM     67  O   TRP A 198      43.482   9.071  11.488  1.00  7.86           O
ATOM     68  CB  TRP A 198      42.033  11.221  13.131  1.00  7.86           C
ATOM     69  CG  TRP A 198      42.767  12.505  12.833  1.00  8.01           C
ATOM     70  CD1 TRP A 198      42.243  13.756  12.803  1.00  8.47           C
ATOM     71  CD2 TRP A 198      44.164  12.653  12.580  1.00  8.09           C
ATOM     72  NE1 TRP A 198      43.219  14.674  12.487  1.00  8.73           N
ATOM     73  CE2 TRP A 198      44.413  14.024  12.350  1.00  8.09           C
ATOM     74  CE3 TRP A 198      45.232  11.768  12.538  1.00  8.16           C
ATOM     75  CZ2 TRP A 198      45.708  14.510  12.050  1.00  8.10           C
ATOM     76  CZ3 TRP A 198      46.496  12.242  12.231  1.00  8.18           C
ATOM     77  CH2 TRP A 198      46.715  13.600  12.019  1.00  8.23           C
ATOM     78  N   ARG A 226      44.174   8.773  20.806  1.00  6.26           N
ATOM     79  CA  ARG A 226      45.520   8.155  20.967  1.00  6.36           C
ATOM     80  C   ARG A 226      46.553   9.119  21.539  1.00  6.53           C
ATOM     81  O   ARG A 226      47.569   8.695  22.043  1.00  6.59           O
ATOM     82  CB  ARG A 226      46.021   7.625  19.619  1.00  6.51           C
ATOM     83  CG  ARG A 226      46.663   8.679  18.715  1.00  6.55           C
ATOM     84  CD  ARG A 226      46.846   8.283  17.267  1.00  6.45           C
ATOM     85  NE  ARG A 226      47.551   9.318  16.549  1.00  6.18           N
ATOM     86  CZ  ARG A 226      47.879   9.277  15.272  1.00  6.50           C
ATOM     87  NH1 ARG A 226      47.480   8.244  14.509  1.00  7.33           N
ATOM     88  NH2 ARG A 226      48.549  10.297  14.764  1.00  6.78           N
ATOM     89  N   GLU A 227      46.321  10.418  21.321  1.00  6.60           N
ATOM     90  CA  GLU A 227      47.312  11.490  21.603  1.00  6.80           C
ATOM     91  C   GLU A 227      46.718  12.848  21.199  1.00  6.74           C
ATOM     92  O   GLU A 227      45.932  12.943  20.249  1.00  6.73           O
ATOM     93  CB  GLU A 227      48.569  11.257  20.731  1.00  7.47           C
ATOM     94  CG  GLU A 227      49.771  12.088  21.155  1.00  8.11           C
ATOM     95  CD  GLU A 227      50.424  12.817  20.036  1.00  8.20           C
ATOM     96  OE1 GLU A 227      50.547  12.305  18.895  1.00  8.84           O
ATOM     97  OE2 GLU A 227      50.804  13.987  20.273  1.00  9.21           O
ATOM     98  N   TYR A 228      47.138  13.905  21.881  1.00  6.66           N
ATOM     99  CA  TYR A 228      46.859  15.283  21.420  1.00  6.71           C
ATOM    100  C   TYR A 228      45.348  15.555  21.380  1.00  6.73           C
ATOM    101  O   TYR A 228      44.626  15.121  22.289  1.00  6.90           O
ATOM    102  CB  TYR A 228      47.592  15.604  20.114  1.00  6.93           C
ATOM    103  CG  TYR A 228      47.971  17.038  19.932  1.00  7.09           C
ATOM    104  CD1 TYR A 228      48.848  17.648  20.817  1.00  6.97           C
ATOM    105  CD2 TYR A 228      47.543  17.742  18.828  1.00  6.98           C
ATOM    106  CE1 TYR A 228      49.272  18.954  20.592  1.00  7.05           C
ATOM    107  CE2 TYR A 228      47.954  19.054  18.588  1.00  7.38           C
ATOM    108  CZ  TYR A 228      48.822  19.627  19.498  1.00  8.06           C
ATOM    109  OH  TYR A 228      49.252  20.952  19.309  1.00  8.46           O
TER
HETATM  110  C1 ANGA B   1      52.184  16.633  15.213  0.50 10.82           C
HETATM  111  C2 ANGA B   1      51.328  15.445  15.648  0.50 10.31           C
HETATM  112  C3 ANGA B   1      50.232  15.929  16.563  0.50  9.54           C
HETATM  113  C4 ANGA B   1      50.840  16.728  17.712  0.50 10.02           C
HETATM  114  C5 ANGA B   1      51.754  17.850  17.194  0.50 10.27           C
HETATM  115  C6 ANGA B   1      52.419  18.596  18.349  0.50 10.78           C
HETATM  116  C7 ANGA B   1      51.064  13.601  14.053  0.50 10.97           C
HETATM  117  C8 ANGA B   1      50.412  13.261  12.752  0.50 11.18           C
HETATM  118  N2 ANGA B   1      50.777  14.837  14.457  0.50 10.28           N
HETATM  119  O1 ANGA B   1      53.223  16.116  14.399  0.50 13.50           O
HETATM  120  O3 ANGA B   1      49.470  14.809  17.068  0.50  7.48           O
HETATM  121  O4 ANGA B   1      51.582  15.870  18.572  0.50 10.21           O
HETATM  122  O5 ANGA B   1      52.747  17.304  16.339  0.50 10.70           O
HETATM  123  O6 ANGA B   1      52.951  19.847  17.881  0.50 11.05           O
HETATM  124  O7 ANGA B   1      51.809  12.830  14.679  0.50 12.42           O
HETATM  125  C1 BA2G B   1      52.458  16.255  15.225  0.50  6.53           C
HETATM  126  C2 BA2G B   1      51.405  15.203  15.602  0.50  6.62           C
HETATM  127  C3 BA2G B   1      50.303  15.813  16.446  0.50  6.74           C
HETATM  128  C4 BA2G B   1      50.898  16.605  17.603  0.50  6.75           C
HETATM  129  C5 BA2G B   1      51.801  17.682  17.009  0.50  6.75           C
HETATM  130  C6 BA2G B   1      52.362  18.579  18.090  0.50  7.01           C
HETATM  131  C7 BA2G B   1      51.061  13.291  14.119  0.50  7.43           C
HETATM  132  C8 BA2G B   1      50.576  12.828  12.782  0.50  7.73           C
HETATM  133  N2 BA2G B   1      50.868  14.598  14.391  0.50  6.93           N
HETATM  134  O1 BA2G B   1      51.913  17.115  14.227  0.50  6.16           O
HETATM  135  O3 BA2G B   1      49.470  14.809  17.068  0.50  7.48           O
HETATM  136  O4 BA2G B   1      51.678  15.787  18.475  0.50  7.10           O
HETATM  137  O5 BA2G B   1      52.868  17.027  16.346  0.50  6.65           O
HETATM  138  O6 BA2G B   1      53.328  19.464  17.495  0.50  7.14           O
HETATM  139  O7 BA2G B   1      51.651  12.510  14.891  0.50  8.28           O
HETATM  140  C1  GAL B   2      48.317  14.425  16.335  1.00  7.03           C
HETATM  141  C2  GAL B   2      47.709  13.209  17.044  1.00  6.63           C
HETATM  142  C3  GAL B   2      46.407  12.813  16.374  1.00  6.30           C
HETATM  143  C4  GAL B   2      45.528  14.023  16.145  1.00  6.50           C
HETATM  144  C5  GAL B   2      46.290  15.197  15.510  1.00  6.90           C
HETATM  145  C6  GAL B   2      45.452  16.426  15.354  1.00  7.16           C
HETATM  146  O2  GAL B   2      48.594  12.104  17.056  1.00  6.77           O
HETATM  147  O3  GAL B   2      45.706  11.865  17.180  1.00  6.72           O
HETATM  148  O4  GAL B   2      44.917  14.359  17.391  1.00  6.48           O
HETATM  149  O5  GAL B   2      47.392  15.490  16.361  1.00  6.90           O
HETATM  150  O6  GAL B   2      46.180  17.431  14.607  1.00  8.08           O
HETATM  151  O   HOH A 402      57.128  11.582  10.607  1.00 37.41           O
HETATM  152  O   HOH A 411      52.843  10.296  14.419  1.00 11.30           O
HETATM  153  O   HOH A 427      12.140   9.626   5.000  1.00 13.06           O
HETATM  154  O   HOH A 428       6.271  10.855  14.160  1.00 22.13           O
HETATM  155  O   HOH A 432      53.837  18.561  13.105  1.00 11.71           O
HETATM  156  O   HOH A 435      52.791  12.317  17.459  1.00 11.07           O
HETATM  157  O   HOH A 437      59.801   5.040  13.918  1.00 18.21           O
HETATM  158  O   HOH A 465      54.179  14.751  18.063  1.00 15.00           O
HETATM  159  O   HOH A 468      51.565  21.790  20.530  1.00 15.45           O
HETATM  160  O   HOH A 492      55.498  18.792  15.849  1.00 36.08           O
HETATM  161  O   HOH A 497      49.326   9.494  12.172  1.00  9.02           O
HETATM  162  O   HOH A 501       5.000   5.247  11.295  1.00 12.56           O
HETATM  163  O   HOH A 534      53.660  20.077  21.205  1.00 20.82           O
HETATM  164  O  BHOH A 553      10.800  11.452   8.547  0.50 21.63           O
HETATM  165  O   HOH A 561      55.104  10.811  17.827  1.00 13.93           O
HETATM  166  O  AHOH A 566      12.525  11.257   7.358  0.50 19.22           O
HETATM  167  O   HOH A 594      55.378  10.746  14.920  1.00 25.54           O
HETATM  168  O  AHOH A 596      56.056  15.010  15.909  0.50 22.86           O
HETATM  169  O  BHOH A 596      55.819  16.141  16.535  0.50 25.09           O
HETATM  170  O   HOH A 604      55.413  15.482  20.526  1.00 20.46           O
HETATM  171  O   HOH A 612      55.669  18.124  20.272  1.00 30.96           O
HETATM  172  O   HOH A 628      52.150   8.845  12.227  0.50  9.11           O
END
'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
