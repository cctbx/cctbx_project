from __future__ import division
from __future__ import print_function
from builtins import range
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
#import libtbx.load_env
from cctbx import geometry_restraints
import iotbx

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

raw_records2 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
HETATM 1406  O   HOH A 282      32.366  19.942  24.727  1.00 38.09           O
END
""".splitlines()

# connect to O, distance 2.9A with symop -y+1,x-y,z-1/3
raw_records3 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
HETATM 1406  O   HOH A 282      32.366  19.942  24.727  1.00 38.09           O
HETATM 1407  O   HOH A 283      33.366  18.942  23.727  1.00 38.09           O
END
""".splitlines()

raw_records4 = """\
CRYST1   15.775   12.565   13.187  90.00  90.00  90.00 P 1
ATOM      1  N   MET A   1       9.821   6.568   5.000  1.00 66.07           N
ATOM      2  CA  MET A   1       9.946   7.171   6.357  1.00 66.55           C
ATOM      3  C   MET A   1      10.571   6.157   7.305  1.00 64.57           C
ATOM      4  O   MET A   1      10.775   5.000   6.933  1.00 66.25           O
ATOM      5  CB  MET A   1       8.570   7.565   6.897  1.00 69.08           C
ATOM      6  CG  MET A   1       7.723   6.373   7.299  1.00 71.37           C
ATOM      7  SD  MET A   1       6.247   6.862   8.187  1.00 76.22           S
ATOM      8  CE  MET A   1       5.000   6.694   6.892  1.00 74.93           C
END
""".splitlines()

raw_records5 = """\
CRYST1  258.687  258.687   47.103  90.00  90.00 120.00 P 63          6
ATOM    213  N   ILE A  78      87.236 -55.209   0.578  1.00179.51           N
ATOM    321  O   LYS A  93      81.801 -49.470  26.164  1.00197.87           O
"""

raw_records6 = """\
CRYST1  101.940  101.370  203.540  90.00  90.00  90.00 C 1 2 1
ATOM      1  N   GLU F 440      51.717  35.260  96.810  1.00173.09           N
ATOM      2  CA  GLU F 440      50.591  34.319  97.027  1.00158.07           C
ATOM      3  C   GLU F 440      49.903  34.410  98.438  1.00181.12           C
ATOM      4  O   GLU F 440      48.809  34.940  98.522  1.00178.50           O
ATOM      5  N   ILE F 441      50.585  34.020  99.528  1.00197.83           N
ATOM      6  CA  ILE F 441      51.029  34.946 100.496  1.00205.24           C
ATOM      7  C   ILE F 441      50.630  36.474 100.471  1.00215.52           C
ATOM      8  O   ILE F 441      50.805  37.162 101.476  1.00228.70           O
ATOM      9  N   GLY F 442      50.138  37.042  99.430  1.00196.69           N
ATOM     10  CA  GLY F 442      49.342  38.214  98.975  1.00156.00           C
ATOM     11  C   GLY F 442      48.218  38.275  97.877  1.00171.08           C
ATOM     12  O   GLY F 442      47.220  38.941  98.131  1.00175.64           O
ATOM     13  N   ILE F 443      48.268  37.635  96.688  1.00184.50           N
ATOM     14  CA  ILE F 443      47.069  37.544  95.818  1.00185.10           C
ATOM     15  C   ILE F 443      46.104  36.600  96.497  1.00180.57           C
ATOM     16  O   ILE F 443      44.887  36.627  96.262  1.00172.35           O
ATOM     17  N   LEU F 444      46.637  35.791  97.397  1.00186.46           N
ATOM     18  CA  LEU F 444      45.893  35.148  98.444  1.00184.47           C
ATOM     19  C   LEU F 444      45.942  35.880  99.806  1.00195.40           C
ATOM     20  O   LEU F 444      44.848  35.788 100.444  1.00193.90           O
"""
raw_records7 = """\
CRYST1  101.940  101.370  203.540  90.00  90.00  90.00 C 1 2 1
ATOM      1  N   ILE F 437      55.794  33.676  96.336  1.00160.31           N
ATOM      2  CA  ILE F 437      54.818  33.116  97.305  1.00161.41           C
ATOM      3  C   ILE F 437      54.482  34.092  98.461  1.00184.09           C
ATOM      4  O   ILE F 437      53.584  33.767  99.360  1.00201.20           O
ATOM      5  N   ASN F 438      55.255  35.246  98.488  1.00184.69           N
ATOM      6  CA  ASN F 438      55.094  36.686  98.818  1.00185.32           C
ATOM      7  C   ASN F 438      53.768  37.307  98.414  1.00190.11           C
ATOM      8  O   ASN F 438      52.882  37.698  99.302  1.00206.13           O
ATOM      9  N   SER F 439      53.683  37.260  97.042  1.00171.17           N
ATOM     10  CA  SER F 439      52.551  37.517  96.151  1.00180.82           C
ATOM     11  C   SER F 439      51.492  36.418  96.171  1.00176.83           C
ATOM     12  O   SER F 439      50.456  36.620  95.555  1.00163.02           O
ATOM     13  N   GLU F 440      51.717  35.260  96.810  1.00173.09           N
ATOM     14  CA  GLU F 440      50.591  34.319  97.027  1.00158.07           C
ATOM     15  C   GLU F 440      49.903  34.410  98.438  1.00181.12           C
ATOM     16  O   GLU F 440      48.809  34.940  98.522  1.00178.50           O
ATOM     17  N   ILE F 441      50.585  34.020  99.528  1.00197.83           N
ATOM     18  CA  ILE F 441      51.029  34.946 100.496  1.00205.24           C
ATOM     19  C   ILE F 441      50.630  36.474 100.471  1.00215.52           C
ATOM     20  O   ILE F 441      50.805  37.162 101.476  1.00228.70           O
"""
# excluded:
#ATOM      1  N   MET A   1       9.821   6.568   5.000  1.00 66.07           N

def make_initial_grm(mon_lib_srv, ener_lib, records):
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = records,
    force_symmetry = True)

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = True,
    plain_pairs_radius = 5.0)
  xrs = processed_pdb_file.xray_structure()
  return geometry, xrs

def show_sorted_geometry(geometry, xrs, file_name):
  out = open(file_name,'w')
  geometry.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=out)
  out.close()

def exercise_add_new_bond_restraint_in_place(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records4)

  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(0,3),
    distance_ideal=2.0,
    weight=3000)
  assert not geometry.is_bonded_atoms(0,3)
  assert not geometry.is_bonded_atoms(3,0)
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart())
  assert geometry.is_bonded_atoms(0,3)
  assert geometry.is_bonded_atoms(3,0)
  assert geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  # That's the way to get them:
  simple, asu = geometry.get_covalent_bond_proxies()
  assert simple.size() + asu.size() == 8
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 8
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0

def exercise_add_super_long_bond(mon_lib_srv, ener_lib):
  # distance between the two is 26A, they are not added because of
  # max_distance_between_connecting_atoms=5 is default.
  # Inspired by 4c8q, atoms are from offending SHEET with wrong parallel/
  # antiparallel definition.
  #
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records5)
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(0,1),
    distance_ideal=2.0,
    weight=3000)
  assert not geometry.is_bonded_atoms(0,1)
  assert not geometry.is_bonded_atoms(1,0)
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
      max_distance_between_connecting_atoms=5)
  assert not geometry.is_bonded_atoms(0,1)
  assert not geometry.is_bonded_atoms(1,0)
  # !!! This will fail, but should not. Left for future investigation.
  # geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
  #     max_distance_between_connecting_atoms=30)

def exercise_bond_near_symmetry(mon_lib_srv, ener_lib):
  """ Making bond near symmetry mate:
  " N   LEU F 444 " -- " O   GLU F 440 " """
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records6)
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(16,3),
    distance_ideal=2.9,
    weight=3000)
  assert not geometry.is_bonded_atoms(16,3)
  assert not geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'before_xercise_bond_near_symmetry.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 19
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
      max_distance_between_connecting_atoms=10)
  assert geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'after_xercise_bond_near_symmetry.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 20  # assert 20
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0  # assert 0

def exercise_bond_near_symmetry2(mon_lib_srv, ener_lib):
  """ Same as previous, other atoms were still failing
  Making bond near symmetry mate:
  " N   LEU F 441 " -- " O   GLU F 437 " """
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records7)
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(16,3),
    distance_ideal=2.9,
    weight=3000)
  assert not geometry.is_bonded_atoms(16,3)
  assert not geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'before_xercise_bond_near_symmetry2.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 19
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
      max_distance_between_connecting_atoms=10)
  assert geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'after_xercise_bond_near_symmetry2.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 20  # assert 20
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0  # assert 0

def exercise_single_atom(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)

  # output for debugging!!!
  # show_sorted_geometry(geometry, xrs, 'before_exersice_single_atoms.geo')

  xrs_add = iotbx.pdb.input(source_info=None, lines=raw_records2) \
      .xray_structure_simple()
  proxy1 = geometry_restraints.bond_simple_proxy(
      i_seqs=(3,9),
      distance_ideal=2.0,
      weight=3000,
      origin_id=1)

  new_xrs = xrs.concatenate(xrs_add)
  all_sites_cart = new_xrs.sites_cart()
  number_of_new_atoms = len(xrs_add.sites_cart())
  new_geometry = geometry.new_included_bonded_atoms(
      proxies=[proxy1],
      sites_cart=all_sites_cart,
      site_symmetry_table=xrs_add.site_symmetry_table(),
      nonbonded_types=flex.std_string(["OH2"]*number_of_new_atoms),
      nonbonded_charges=flex.int(number_of_new_atoms, 0),
      skip_max_proxy_distance_calculation=True)
  # output for debugging!!!
  # show_sorted_geometry(new_geometry, new_xrs,
  #     'after_exersice_single_atoms.geo')
  assert new_geometry.is_bonded_atoms(3,9)

  assert new_geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert new_geometry.pair_proxies().bond_proxies.asu.size() == 1
  # That's the way to get them:
  simple, asu = new_geometry.get_covalent_bond_proxies()
  assert simple.size() + asu.size() == 8
  simple, asu = new_geometry.get_all_bond_proxies()
  assert simple.size() + asu.size() == 9, "%d, %d" % (simple.size(), asu.size())
  assert new_geometry.pair_proxies().nonbonded_proxies.simple.size() == 10
  assert new_geometry.pair_proxies().nonbonded_proxies.asu.size() ==2
  assert new_geometry.get_hbond_proxies_iseqs() == [(3, 9)]
  simple, asu = new_geometry.get_covalent_bond_proxies()
  assert simple.size() + asu.size() == 8, "%d" % (simple.size() + asu.size())

def exercise_multiple_atoms(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)

  # output for debugging!!!
  # show_sorted_geometry(geometry, xrs, 'before_exersice_multiple_atoms.geo')

  xrs_add = iotbx.pdb.input(source_info=None, lines=raw_records3) \
      .xray_structure_simple()
  proxy1 = geometry_restraints.bond_simple_proxy(
      i_seqs=(3,9),
      distance_ideal=2.0,
      weight=3000)
  proxy2 = geometry_restraints.bond_simple_proxy(
      i_seqs=(4,10),
      distance_ideal=2.0,
      weight=3000)

  new_xrs = xrs.concatenate(xrs_add)
  all_sites_cart = new_xrs.sites_cart()
  number_of_new_atoms = len(xrs_add.sites_cart())
  assert geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 10
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0
  new_geometry = geometry.new_included_bonded_atoms(
      proxies=[proxy1, proxy2],
      sites_cart=all_sites_cart,
      site_symmetry_table=xrs_add.site_symmetry_table(),
      nonbonded_types=flex.std_string(["OH2"]*number_of_new_atoms),
      nonbonded_charges=flex.int(number_of_new_atoms, 0),
      skip_max_proxy_distance_calculation=True)
  # output for debugging!!!
  # show_sorted_geometry(new_geometry, new_xrs,
  #     'after_exersice_multiple_atoms.geo')
  assert new_geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert new_geometry.pair_proxies().bond_proxies.asu.size() == 2
  assert new_geometry.pair_proxies().nonbonded_proxies.simple.size() == 11
  assert new_geometry.pair_proxies().nonbonded_proxies.asu.size() == 4

def exercise_is_bonded_atoms(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)
  # show_sorted_geometry(geometry, xrs, 'blabla.geo')
  linked_atoms = [(0,1), (1,2), (1,4), (2,3), (4,5),(5,6), (6,7), (7,8)]
  for i in range (9):
    for j in range(9):
      assert geometry.is_bonded_atoms(i,j) == geometry.is_bonded_atoms(j,i)
      assert geometry.is_bonded_atoms(i,j) == (tuple(sorted([i,j])) in linked_atoms)

def exercise():
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None and ener_lib is not None:
    exercise_add_new_bond_restraint_in_place(mon_lib_srv, ener_lib)
    exercise_single_atom(mon_lib_srv, ener_lib)
    exercise_multiple_atoms(mon_lib_srv, ener_lib)
    exercise_is_bonded_atoms(mon_lib_srv, ener_lib)
    exercise_add_super_long_bond(mon_lib_srv, ener_lib)
    exercise_bond_near_symmetry(mon_lib_srv, ener_lib)
    exercise_bond_near_symmetry2(mon_lib_srv, ener_lib)

if (__name__ == "__main__"):
  exercise()
