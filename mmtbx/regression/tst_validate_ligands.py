from __future__ import absolute_import, division, print_function
import time, os
from libtbx.utils import null_out
import libtbx.load_env
from libtbx.test_utils import approx_equal
from iotbx.cli_parser import run_program
from mmtbx.programs import validate_ligands as val_lig

# ------------------------------------------------------------------------------

def run():
  run_test1()
  run_test2()
  run_test3()
  run_test4()

# ------------------------------------------------------------------------------

def run_test1():
  '''
  Test if iselection for ligand PG5 (chain A resseq 201) is correct.
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
  for lr in vl_manager:
    assert (lr.id_str == 'PG5 A 201')
  # test iselection
  assert(list(lr.ligand_isel) == [190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
    200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214,
    215, 216, 217, 218, 219])

  for lr in vl_manager:
    clashes_result = lr.get_overlaps()
    assert(clashes_result.n_clashes == 5)
    #print(round(clashes_result.clashscore,1), 27.6)
    assert approx_equal(clashes_result.clashscore, 27.6, eps=0.5)

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

  os.remove('two_chains_ligand_water_newH.cif')
  os.remove('two_chains_ligand_water_newH.txt')

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
      #print(round(clashes_result.clashscore,1))
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

      assert approx_equal(adps.b_min_within, 4.00, eps=0.01)
      assert approx_equal(adps.b_max_within, 54.65, eps=0.01)
      assert approx_equal(adps.b_mean_within, 23.23, eps=0.02)
      #
      assert(clashes_result.n_clashes == 4)
      #print(round(clashes_result.clashscore,1))
      assert approx_equal(clashes_result.clashscore, 13.0, eps=0.5)
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
      #print(round(clashes_result.clashscore,1))
      assert approx_equal(clashes_result.clashscore, 18.5, eps=0.5)

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

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
