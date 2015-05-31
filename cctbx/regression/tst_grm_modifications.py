from __future__ import division
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
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 8
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0

def exercise_single_atom(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)

  # output for debugging!!!
  # show_sorted_geometry(geometry, xrs, 'before_exersice_single_atoms.geo')

  xrs_add = iotbx.pdb.input(source_info=None, lines=raw_records2) \
      .xray_structure_simple()
  proxy1 = geometry_restraints.bond_simple_proxy(
      i_seqs=(3,9),
      distance_ideal=2.0,
      weight=3000)

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
  assert new_geometry.pair_proxies().nonbonded_proxies.simple.size() == 10
  assert new_geometry.pair_proxies().nonbonded_proxies.asu.size() ==2

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
  show_sorted_geometry(geometry, xrs, 'blabla.geo')
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
    print "Can not initialize monomer_library, skipping test."
  if mon_lib_srv is not None and ener_lib is not None:
    exercise_add_new_bond_restraint_in_place(mon_lib_srv, ener_lib)
    exercise_single_atom(mon_lib_srv, ener_lib)
    exercise_multiple_atoms(mon_lib_srv, ener_lib)
    exercise_is_bonded_atoms(mon_lib_srv, ener_lib)

if (__name__ == "__main__"):
  exercise()
