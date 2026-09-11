from __future__ import absolute_import, division, print_function
import iotbx.pdb, iotbx.cif
from mmtbx.monomer_library import server
from mmtbx.hydrogens import reduce_hydrogen

def run():
  test_001()
  test_002()
  test_003()

# ------------------------------------------------------------------------------

def test_001():
  '''
    A restraint file with several ligands: each ligand's CH2 names are checked
    against its own ideal coordinates. AC3 and AD0 (codes not in the CCD) share
    atom names, with H21/H22 swapped between them. The AD0 H already match AD0,
    so name_prochiral_h must leave them; it used to read AC3, the first block.
  '''
  file_name = "tst_add_hydrogen_16.cif"
  with open(file_name, "w") as f:
    f.write(cif_str_001)
  mon_lib_srv = server.server()
  mon_lib_srv.process_cif_object(
    cif_object = iotbx.cif.reader(input_string=cif_str_001).model(),
    file_name  = file_name)
  h = iotbx.pdb.input(lines=pdb_str_001.split("\n"), source_info=None).construct_hierarchy()
  h.atoms().reset_i_seq()
  names = list(h.atoms().extract_name())
  reduce_hydrogen.name_prochiral_h(h, mon_lib_srv)
  assert list(h.atoms().extract_name()) == names, (
    names, list(h.atoms().extract_name()))

def get_user_lig_srv():
  '''User restraint file for LIG, a code the CCD uses for another molecule.'''
  file_name = "tst_add_hydrogen_16_lig.cif"
  with open(file_name, "w") as f:
    f.write(cif_str_002)
  mon_lib_srv = server.server()
  mon_lib_srv.process_cif_object(
    cif_object = iotbx.cif.reader(input_string=cif_str_002).model(),
    file_name  = file_name)
  return mon_lib_srv

def test_002():
  '''
    A user restraint file wins over an unrelated CCD entry with the same code.
    CCD LIG (C15H11N3) has a CH2 at C17 with the same names as this user LIG
    (propane C5-C17-C20) but the opposite handedness; the H match the user file,
    so name_prochiral_h must leave them (it used CCD geometry and swapped them).
  '''
  mon_lib_srv = get_user_lig_srv()
  h = iotbx.pdb.input(lines=pdb_str_002.split("\n"), source_info=None).construct_hierarchy()
  h.atoms().reset_i_seq()
  names = list(h.atoms().extract_name())
  reduce_hydrogen.name_prochiral_h(h, mon_lib_srv)
  assert list(h.atoms().extract_name()) == names, (
    names, list(h.atoms().extract_name()))

def test_003():
  '''
    Bond orders come from the user file only: the user LIG has single bonds, CCD
    LIG's double bonds (C5=C6, C20=C21, ...) must not leak in.
  '''
  mon_lib_srv = get_user_lig_srv()
  orders = reduce_hydrogen._bond_orders("LIG", mon_lib_srv, {})
  assert orders == {}, sorted(tuple(sorted(k)) for k in orders)

# ------------------------------------------------------------------------------

cif_str_001 = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
AC3 AC3 'ethanol-like, first ligand' ligand 9 3 .
AD0 AD0 'ethanol-like, second ligand' ligand 9 3 .

data_comp_AC3
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
   AC3 C1   C   1.5200   0.0000   0.0000
   AC3 C2   C   0.0000   0.0000   0.0000
   AC3 O3   O  -0.4773   1.3480   0.0000
   AC3 H11  H   1.8641   0.8604   0.2868
   AC3 H12  H   1.8678  -0.5797   0.6956
   AC3 H13  H   1.8641  -0.2868  -0.8604
   AC3 H21  H  -0.3231  -0.4572  -0.7921
   AC3 H22  H  -0.3231  -0.4572   0.7921
   AC3 H3   H  -1.2837   1.6168   0.0000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
   AC3 C1   C2   single 1.520 0.020
   AC3 C2   O3   single 1.430 0.020
   AC3 C1   H11  single 0.970 0.020
   AC3 C1   H12  single 0.970 0.020
   AC3 C1   H13  single 0.970 0.020
   AC3 C2   H21  single 0.970 0.020
   AC3 C2   H22  single 0.970 0.020
   AC3 O3   H3   single 0.850 0.020

data_comp_AD0
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
   AD0 C1   C   1.5200   0.0000   0.0000
   AD0 C2   C   0.0000   0.0000   0.0000
   AD0 O3   O  -0.4773   1.3480   0.0000
   AD0 H11  H   1.8641   0.8604   0.2868
   AD0 H12  H   1.8678  -0.5797   0.6956
   AD0 H13  H   1.8641  -0.2868  -0.8604
   AD0 H21  H  -0.3231  -0.4572   0.7921
   AD0 H22  H  -0.3231  -0.4572  -0.7921
   AD0 H3   H  -1.2837   1.6168   0.0000
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
   AD0 C1   C2   single 1.520 0.020
   AD0 C2   O3   single 1.430 0.020
   AD0 C1   H11  single 0.970 0.020
   AD0 C1   H12  single 0.970 0.020
   AD0 C1   H13  single 0.970 0.020
   AD0 C2   H21  single 0.970 0.020
   AD0 C2   H22  single 0.970 0.020
   AD0 O3   H3   single 0.850 0.020
"""

pdb_str_001 = """
HETATM    1  C1  AD0 A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    2  C2  AD0 A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    3  O3  AD0 A   1      -0.477   1.348   0.000  1.00 20.00           O
HETATM    4  H21 AD0 A   1      -0.323  -0.457   0.792  1.00 20.00           H
HETATM    5  H22 AD0 A   1      -0.323  -0.457  -0.792  1.00 20.00           H
END
"""

cif_str_002 = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LIG LIG 'user ligand, propane' ligand 11 3 .

data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
   LIG C5   C  -0.5694   1.4093   0.0000
   LIG C17  C   0.0000   0.0000   0.0000
   LIG C20  C   1.5200   0.0000   0.0000
   LIG H171 H  -0.3131  -0.4641   0.7921
   LIG H172 H  -0.3131  -0.4641  -0.7921
   LIG H5A  H   0.1571   2.0521   0.0000
   LIG H5B  H  -1.1146   1.5383  -0.7919
   LIG H5C  H  -1.1146   1.5383   0.7919
   LIG H20A H   1.8438  -0.9144   0.0000
   LIG H20B H   1.8438   0.4572  -0.7919
   LIG H20C H   1.8438   0.4572   0.7919
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
   LIG C5   C17  single 1.520 0.020
   LIG C17  C20  single 1.520 0.020
   LIG C17  H171 single 0.970 0.020
   LIG C17  H172 single 0.970 0.020
   LIG C5   H5A  single 0.970 0.020
   LIG C5   H5B  single 0.970 0.020
   LIG C5   H5C  single 0.970 0.020
   LIG C20  H20A single 0.970 0.020
   LIG C20  H20B single 0.970 0.020
   LIG C20  H20C single 0.970 0.020
"""

pdb_str_002 = """
HETATM    1  C5  LIG A   1      -0.569   1.409   0.000  1.00 20.00           C
HETATM    2  C17 LIG A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    3  C20 LIG A   1       1.520   0.000   0.000  1.00 20.00           C
HETATM    4 H171 LIG A   1      -0.313  -0.464   0.792  1.00 20.00           H
HETATM    5 H172 LIG A   1      -0.313  -0.464  -0.792  1.00 20.00           H
END
"""

if __name__ == '__main__':
  run()
  print("OK")
