from __future__ import absolute_import, division, print_function
import iotbx.pdb, iotbx.cif
import mmtbx.model
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from mmtbx.monomer_library import server
from mmtbx.hydrogens import reduce_hydrogen

def run():
  test_001()
  test_002()
  test_003()
  test_004()
  test_005()
  test_006()
  test_007()
  test_008()

# ------------------------------------------------------------------------------

def place(pdb_str):
  '''Place H (no optimization); return the atoms by name.'''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  return {a.name.strip(): a for a in obj.get_model().get_hierarchy().atoms()}

def h_names_on(pdb_str, parent_name):
  '''Place H (no optimization); return names of the H bonded to parent_name.'''
  atoms = list(place(pdb_str).values())
  parent = [a for a in atoms if a.name.strip() == parent_name][0]
  return sorted(a.name.strip() for a in atoms
    if a.element.strip() == 'H' and a.distance(parent) < 1.2)

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

def test_004():
  '''
    6EL (5yj1) passes test_for_peptide (N, CA, C, O) but its H2 sits on ring
    CG, not on N. The N-terminal H2 filter removed it by name: CG got 1 H.
  '''
  names = h_names_on(pdb_str_004, 'CG')
  assert names == ['H2', 'H3'], names

def test_005():
  '''
    Control: 0TD H2 is the amine H on N; a residue that is not first in the
    chain keeps only H on N.
  '''
  names = h_names_on(pdb_str_005, 'N')
  assert names == ['H'], names

def test_006():
  '''
    GLZ (aminoacetaldehyde) passes test_for_peptide; its HXT is the aldehyde
    H on C (no OXT). The C-terminal HXT filter removed it by name.
  '''
  names = h_names_on(pdb_str_006, 'C')
  assert names == ['HXT'], names

def check_h(atoms, h, a0, a1, angle_ideal):
  '''H sits on a0 at a plausible X-H length and H-a0-a1 angle.'''
  d = atoms[h].distance(atoms[a0])
  assert 0.8 < d < 1.15, (h, d)
  angle = atoms[a0].angle(atoms[h], atoms[a1], deg=True)
  assert approx_equal(angle, angle_ideal, eps=1.0), (h, angle)

def test_007():
  '''
    PEO (H2O2): nothing beyond the O-O bond anchors the H dihedral, so riding H
    cannot parameterize HO1/HO2 and they were deleted (1ng4).
  '''
  atoms = place(pdb_str_007)
  assert sorted(n for n in atoms if n.startswith('H')) == ['HO1', 'HO2'], list(atoms)
  check_h(atoms, 'HO1', 'O1', 'O2', 100.92)
  check_h(atoms, 'HO2', 'O2', 'O1', 100.89)

def test_008():
  '''
    MOH (methanol): same for the methyl and OH H. EOH (ethanol, anchored by the
    third heavy atom) is the control.
  '''
  atoms = place(pdb_str_008)
  hs = sorted(n for n in atoms if n.startswith('H'))
  assert hs == ['H1', 'H2', 'H3', 'HO'], hs
  for h, angle_ideal in [('H1', 107.79), ('H2', 112.98), ('H3', 112.94)]:
    check_h(atoms, h, 'C', 'O', angle_ideal)
  check_h(atoms, 'HO', 'O', 'C', 108.03)
  for h1, h2 in [('H1', 'H2'), ('H1', 'H3'), ('H2', 'H3')]:
    assert atoms[h1].distance(atoms[h2]) > 1.4, (h1, h2)
  hs = sorted(n for n in place(pdb_str_008_control) if n.startswith('H'))
  assert hs == ['H11', 'H12', 'H21', 'H22', 'H23', 'HO'], hs

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

pdb_str_004 = """
CRYST1  201.615  201.615  123.572  90.00  90.00 120.00 H 3
HETATM    1  N   6EL A 501     -18.527  39.011   4.685  1.00 42.85           N
HETATM    2  CA  6EL A 501     -19.238  37.942   3.975  1.00 41.65           C
HETATM    3  C   6EL A 501     -19.091  38.062   2.480  1.00 39.05           C
HETATM    4  O   6EL A 501     -18.732  39.122   1.981  1.00 35.88           O
HETATM    5  CB  6EL A 501     -18.971  36.528   4.429  1.00 41.97           C
HETATM    6  CG  6EL A 501     -20.045  35.653   3.781  1.00 34.50           C
HETATM    7  CD  6EL A 501     -20.198  35.965   2.286  1.00 34.24           C
HETATM    8  OE1 6EL A 501     -20.787  35.166   1.567  1.00 32.94           O
HETATM    9  NE2 6EL A 501     -19.708  37.135   1.735  1.00 37.42           N
HETATM   10  CAE 6EL A 501     -16.110  42.317   6.704  1.00 46.24           C
HETATM   11  CAF 6EL A 501     -17.290  43.026   6.465  1.00 43.30           C
HETATM   12  CAG 6EL A 501     -15.983  40.985   6.304  1.00 44.16           C
HETATM   13  CAH 6EL A 501     -18.362  42.399   5.839  1.00 46.24           C
HETATM   14  CAN 6EL A 501     -17.263  39.108   5.113  1.00 43.31           C
HETATM   15  CAO 6EL A 501     -19.045  40.227   4.740  1.00 46.27           C
HETATM   16  CAP 6EL A 501     -17.059  40.369   5.678  1.00 43.68           C
HETATM   17  CAQ 6EL A 501     -18.247  41.059   5.505  1.00 41.14           C
HETATM   18  OAC 6EL A 501     -16.455  38.154   5.201  1.00 37.79           O
HETATM   19  OAD 6EL A 501     -20.260  40.446   4.509  1.00 36.94           O
END
"""

pdb_str_005 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   0TD A 501       9.131  12.228   9.576  1.00 20.00           N
HETATM    2  CA  0TD A 501       9.235  10.781   9.451  1.00 20.00           C
HETATM    3  C   0TD A 501       7.932  10.147  10.006  1.00 20.00           C
HETATM    4  O   0TD A 501       7.556  10.575  11.117  1.00 20.00           O
HETATM    5  CSB 0TD A 501      13.190  10.241  10.653  1.00 20.00           C
HETATM    6  SB  0TD A 501      11.968  10.858   9.476  1.00 20.00           S
HETATM    7  CB  0TD A 501      10.441  10.157  10.187  1.00 20.00           C
HETATM    8  CG  0TD A 501      10.431   8.611  10.102  1.00 20.00           C
HETATM    9  OD2 0TD A 501      10.667   8.104   8.990  1.00 20.00           O
HETATM   10  OD1 0TD A 501      10.165   8.013  11.163  1.00 20.00           O
HETATM   11  OXT 0TD A 501       7.383   9.269   9.316  1.00 20.00           O
END
"""

pdb_str_006 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   GLZ A 501       9.839   9.850  11.172  1.00 20.00           N
HETATM    2  CA  GLZ A 501       8.601   9.331  10.629  1.00 20.00           C
HETATM    3  C   GLZ A 501       8.696   7.841  10.440  1.00 20.00           C
HETATM    4  O   GLZ A 501       8.048   7.204   9.648  1.00 20.00           O
END
"""

pdb_str_007 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  O1  PEO A 501      14.917  15.127  15.671  1.00 20.00           O
HETATM    2  O2  PEO A 501      15.083  15.207  14.264  1.00 20.00           O
END
"""

pdb_str_008 = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  C   MOH A 501      12.578  10.567  10.234  1.00 20.00           C
HETATM    2  O   MOH A 501      13.940  10.360   9.945  1.00 20.00           O
END
"""

pdb_str_008_control = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EOH A 501      10.157  14.225  12.482  1.00 20.00           C
HETATM    2  C2  EOH A 501      11.073  15.376  12.125  1.00 20.00           C
HETATM    3  O   EOH A 501      10.877  13.288  13.257  1.00 20.00           O
END
"""

if __name__ == '__main__':
  run()
  print("OK")
