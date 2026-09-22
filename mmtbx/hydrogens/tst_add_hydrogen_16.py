from __future__ import absolute_import, division, print_function
import iotbx.pdb, iotbx.cif
import mmtbx.model
from scitbx import matrix
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
  test_009()
  test_010()

# ------------------------------------------------------------------------------

def place_atoms(pdb_str):
  '''Place H (no optimization); return all atoms. Names repeat between
  residues, so only single-residue models may be keyed by name (place()).'''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  return list(obj.get_model().get_hierarchy().atoms())

def place(pdb_str):
  '''Place H in a one-residue model; return the atoms by name.'''
  return {a.name.strip(): a for a in place_atoms(pdb_str)}

def h_names_on(pdb_str, parent_name):
  '''Place H (no optimization); return names of the H bonded to parent_name.'''
  atoms = place_atoms(pdb_str)
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

def h_on(atoms, resname, resseq, parent_name):
  '''Names of the H bonded to <parent_name> of one residue.'''
  parent = [a for a in atoms
    if a.name.strip() == parent_name
    and a.parent().resname.strip() == resname
    and a.parent().parent().resseq.strip() == resseq]
  assert len(parent) == 1, (resname, resseq, parent_name)
  return sorted(a.name.strip() for a in atoms
    if a.element.strip() == 'H' and a.distance(parent[0]) < 1.3)

def test_009():
  '''
    5nxq: 9FZ and 9G2 are in the polypeptide chain, but their dictionaries are
    the free amino acid (N is NH2: H3/H4 and H6/H7). Both H stayed, so N had
    four bonds. The terminal-H filter only knew the name H2.
  '''
  atoms = place_atoms(pdb_str_009)
  assert len(h_on(atoms, '9FZ', '4', 'N')) == 1, h_on(atoms, '9FZ', '4', 'N')
  assert len(h_on(atoms, '9G2', '10', 'N')) == 1, h_on(atoms, '9G2', '10', 'N')
  # the amino acids in the same chain are unchanged
  assert h_on(atoms, 'ILE', '5', 'N') == ['H'], h_on(atoms, 'ILE', '5', 'N')
  assert h_on(atoms, 'GLU', '11', 'N') == ['H'], h_on(atoms, 'GLU', '11', 'N')
  # the C9-C10 cross-link between the two ligands still removes one methyl H
  assert len(h_on(atoms, '9FZ', '4', 'C9')) == 2, h_on(atoms, '9FZ', '4', 'C9')
  assert len(h_on(atoms, '9G2', '10', 'C10')) == 2, h_on(atoms, '9G2', '10', 'C10')

def ch2_hand(atoms, resname, resseq, parent, hv1, hv2, h2):
  '''
  Sign of the chiral volume (hv1, hv2, h2) at parent: which side of the two
  heavy neighbours the '2' hydrogen of a CH2 is named on.
  '''
  a = dict((x.name.strip(), matrix.col(x.xyz)) for x in atoms
    if x.parent().resname.strip() == resname
    and x.parent().parent().resseq.strip() == resseq)
  p = a[parent]
  return 1 if ((a[hv1]-p).cross(a[hv2]-p)).dot(a[h2]-p) > 0 else -1

def test_010():
  '''
    LEU HB2/HB3 came out swapped (7IT7 fragment). The CH2 reference read the
    CCD's ideal coordinates, and for ARG CB/CG, ILE CG1, LEU CB and MET CB/CG
    those contradict the same entry's model coordinates. The model coordinates
    are what reduce and deposited models follow; all three groups here are -1.
    GLU CB/CG (the CCD agrees with itself) are the control.
  '''
  atoms = place_atoms(pdb_str_010)
  assert ch2_hand(atoms, 'LEU', '107', 'CB', 'CA', 'CG', 'HB2') == -1
  assert ch2_hand(atoms, 'GLU', '41', 'CB', 'CA', 'CG', 'HB2') == -1
  assert ch2_hand(atoms, 'GLU', '41', 'CG', 'CB', 'CD', 'HG2') == -1

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

pdb_str_009 = """
CRYST1   88.679  100.287  219.749  90.00  90.00  90.00 P 2 21 21
ATOM      1  N   ILE D   3      -4.315  36.803  -4.113  1.00 99.90           N
ATOM      2  CA  ILE D   3      -2.949  37.203  -4.414  1.00 90.56           C
ATOM      3  C   ILE D   3      -2.949  38.463  -5.269  1.00 91.77           C
ATOM      4  O   ILE D   3      -3.680  38.556  -6.255  1.00 91.19           O
ATOM      5  CB  ILE D   3      -2.182  36.067  -5.114  1.00 86.55           C
ATOM      6  CG1 ILE D   3      -2.220  34.795  -4.265  1.00 87.21           C
ATOM      7  CG2 ILE D   3      -0.747  36.482  -5.399  1.00 77.11           C
ATOM      8  CD1 ILE D   3      -1.760  33.558  -5.003  1.00 82.40           C
HETATM    9  N   9FZ D   4      -2.130  39.434  -4.879  1.00 91.16           N
HETATM   10  CA  9FZ D   4      -1.988  40.695  -5.655  1.00 89.10           C
HETATM   11  C   9FZ D   4      -1.176  40.466  -6.907  1.00 90.90           C
HETATM   12  O   9FZ D   4       0.010  40.188  -6.831  1.00 91.37           O
HETATM   13  CB  9FZ D   4      -1.337  41.773  -4.798  1.00 89.28           C
HETATM   14  C1  9FZ D   4      -1.144  43.047  -5.607  1.00 90.46           C
HETATM   15  C2  9FZ D   4      -0.679  44.202  -4.730  1.00 97.10           C
HETATM   16  C7  9FZ D   4      -1.020  47.117  -6.816  1.00106.60           C
HETATM   17  C8  9FZ D   4      -1.634  45.994  -6.287  1.00107.15           C
HETATM   18  C9  9FZ D   4      -1.906  47.950  -7.705  1.00102.31           C
HETATM   19  N4  9FZ D   4      -0.689  45.413  -5.545  1.00104.45           N
HETATM   20  N5  9FZ D   4       0.515  46.128  -5.583  1.00105.46           N
HETATM   21  N6  9FZ D   4       0.267  47.232  -6.418  1.00106.68           N
ATOM     22  N   ILE D   5      -1.820  40.576  -8.065  1.00 90.55           N
ATOM     23  CA  ILE D   5      -1.159  40.296  -9.335  1.00 92.65           C
ATOM     24  C   ILE D   5      -1.295  41.456 -10.315  1.00 90.25           C
ATOM     25  O   ILE D   5      -0.639  41.471 -11.356  1.00 88.80           O
ATOM     26  CB  ILE D   5      -1.710  39.006  -9.973  1.00 84.60           C
ATOM     27  CG1 ILE D   5      -3.209  39.150 -10.250  1.00 82.15           C
ATOM     28  CG2 ILE D   5      -1.432  37.800  -9.086  1.00 79.75           C
ATOM     29  CD1 ILE D   5      -3.813  37.968 -10.971  1.00 78.88           C
ATOM     30  N   LEU D   9       0.789  42.265 -13.727  1.00 85.51           N
ATOM     31  CA  LEU D   9       0.169  42.262 -15.045  1.00 87.53           C
ATOM     32  C   LEU D   9       0.491  43.540 -15.814  1.00 93.45           C
ATOM     33  O   LEU D   9       0.677  43.508 -17.030  1.00 95.25           O
ATOM     34  CB  LEU D   9      -1.347  42.092 -14.925  1.00 86.95           C
ATOM     35  CG  LEU D   9      -1.825  40.801 -14.256  1.00 85.95           C
ATOM     36  CD1 LEU D   9      -3.346  40.738 -14.220  1.00 82.56           C
ATOM     37  CD2 LEU D   9      -1.241  39.579 -14.952  1.00 67.49           C
HETATM   38  N   9G2 D  10       0.559  44.662 -15.103  1.00 96.10           N
HETATM   39  CA  9G2 D  10       0.902  45.962 -15.743  1.00 97.17           C
HETATM   40  C   9G2 D  10       2.361  46.013 -16.127  1.00 95.05           C
HETATM   41  O   9G2 D  10       2.735  46.711 -17.055  1.00 95.78           O
HETATM   42  CB  9G2 D  10       0.575  47.130 -14.821  1.00 95.86           C
HETATM   43  C1  9G2 D  10      -0.906  47.477 -14.891  1.00 98.06           C
HETATM   44  C10 9G2 D  10      -1.116  49.129  -8.268  1.00106.65           C
HETATM   45  C2  9G2 D  10      -1.181  48.865 -14.323  1.00107.10           C
HETATM   46  C7  9G2 D  10      -0.869  48.747 -10.699  1.00104.41           C
HETATM   47  C8  9G2 D  10      -0.231  48.723 -11.936  1.00101.61           C
HETATM   48  C9  9G2 D  10      -0.181  48.654  -9.359  1.00104.27           C
HETATM   49  N4  9G2 D  10      -1.186  48.828 -12.861  1.00107.44           N
HETATM   50  N5  9G2 D  10      -2.449  48.916 -12.275  1.00105.64           N
HETATM   51  N6  9G2 D  10      -2.206  48.860 -10.894  1.00106.43           N
ATOM     52  N   GLU D  11       3.190  45.273 -15.398  1.00 91.80           N
ATOM     53  CA  GLU D  11       4.607  45.175 -15.723  1.00 97.96           C
ATOM     54  C   GLU D  11       4.775  44.382 -17.015  1.00 99.71           C
ATOM     55  O   GLU D  11       5.723  44.593 -17.772  1.00103.15           O
ATOM     56  CB  GLU D  11       5.382  44.518 -14.579  1.00104.15           C
ATOM     57  CG  GLU D  11       6.888  44.481 -14.781  1.00112.57           C
ATOM     58  CD  GLU D  11       7.602  43.692 -13.700  1.00128.77           C
ATOM     59  OE1 GLU D  11       6.919  43.144 -12.809  1.00117.73           O
ATOM     60  OE2 GLU D  11       8.848  43.618 -13.744  1.00136.20           O
END
"""

pdb_str_010 = """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A  41      14.715 -15.370  14.923  1.00 39.91           N
ATOM      2  CA  GLU A  41      14.788 -16.780  14.439  1.00 35.65           C
ATOM      3  C   GLU A  41      13.961 -16.899  13.156  1.00 36.25           C
ATOM      4  O   GLU A  41      14.465 -17.527  12.205  1.00 38.69           O
ATOM      5  CB  GLU A  41      14.321 -17.749  15.518  1.00 33.50           C
ATOM      6  CG  GLU A  41      15.333 -17.993  16.600  1.00 35.41           C
ATOM      7  CD  GLU A  41      16.636 -18.604  16.105  1.00 43.11           C
ATOM      8  OE1 GLU A  41      16.630 -19.181  15.035  1.00 45.69           O
ATOM      9  OE2 GLU A  41      17.647 -18.494  16.799  1.00 51.56           O
ATOM     37  N   LEU A 107       4.234 -18.401  11.202  1.00 29.52           N
ATOM     38  CA  LEU A 107       5.506 -19.021  10.766  1.00 31.21           C
ATOM     39  C   LEU A 107       5.985 -20.013  11.824  1.00 29.33           C
ATOM     40  O   LEU A 107       5.845 -19.764  13.033  1.00 32.67           O
ATOM     41  CB  LEU A 107       6.554 -17.940  10.507  1.00 31.28           C
ATOM     42  CG  LEU A 107       6.196 -16.889   9.461  1.00 33.41           C
ATOM     43  CD1 LEU A 107       7.334 -15.910   9.340  1.00 35.68           C
ATOM     44  CD2 LEU A 107       5.884 -17.500   8.088  1.00 32.97           C
END
"""

if __name__ == '__main__':
  run()
  print("OK")
