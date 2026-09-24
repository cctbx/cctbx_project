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
  test_011()
  test_012()
  test_013()

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

def atom_of(atoms, resname, resseq, name):
  '''The one atom of that residue.'''
  found = [a for a in atoms if a.name.strip() == name
    and a.parent().resname.strip() == resname
    and a.parent().parent().resseq.strip() == resseq]
  assert len(found) == 1, (resname, resseq, name)
  return found[0]

def test_013():
  '''
    5nxq again: the amide H of a ligand in the chain. The TRANS link defines
    C-N-H by name, so a residue whose H is called anything else (9FZ H3, 9G2
    H6) got no such restraint, riding called it a rotatable amine H and it
    ended up 17 deg from the preceding C. The amide H now sits where ILE 5's
    does (the control): the three angles at N add up to 360 and it is trans to
    the preceding O. CA-N-H comes from the link as well, since the dictionary
    of the free amino acid gives N the sp3 amine value (9FZ 108.5 deg).
  '''
  atoms = place_atoms(pdb_str_013)
  for resname, resseq, h_name, prev_resname, prev_resseq in (
      ('9FZ', '4',  'H3', 'ILE', '3'),
      ('9G2', '10', 'H6', 'LEU', '9'),
      ('ILE', '5',  'H',  '9FZ', '4')):
    n  = atom_of(atoms, resname, resseq, 'N')
    h  = atom_of(atoms, resname, resseq, h_name)
    ca = atom_of(atoms, resname, resseq, 'CA')
    c  = atom_of(atoms, prev_resname, prev_resseq, 'C')
    o  = atom_of(atoms, prev_resname, prev_resseq, 'O')
    why = (resname, h_name)
    assert approx_equal(n.angle(c, h, deg=True), 124.3, eps=1.0), why
    assert approx_equal(n.angle(ca, h, deg=True), 114.0, eps=1.0), why
    assert approx_equal(n.angle(c, h, deg=True) + n.angle(ca, h, deg=True)
                        + n.angle(c, ca, deg=True), 360.0, eps=1.0), why
    d = matrix.col(h.xyz) - matrix.col(n.xyz)
    plane = ((matrix.col(c.xyz) - matrix.col(n.xyz)).cross(
              matrix.col(ca.xyz) - matrix.col(n.xyz))).normalize()
    assert abs(plane.dot(d)) < 0.1, (why, 'H is out of the amide plane')
    assert 0.8 < n.distance(h) < 1.1, (why, n.distance(h))

def hand(atoms, resname, resseq, parent, first, second, third):
  '''
  Sign of the chiral volume (first, second, third) at parent. For a CH2 pass
  the two heavy neighbours and the '2' hydrogen: which side that H is named on.
  For a propeller pass two H and the heavy neighbour: the turning sense of the
  names around the axis.
  '''
  a = dict((x.name.strip(), matrix.col(x.xyz)) for x in atoms
    if x.parent().resname.strip() == resname
    and x.parent().parent().resseq.strip() == resseq)
  p = a[parent]
  return 1 if ((a[first]-p).cross(a[second]-p)).dot(a[third]-p) > 0 else -1

def test_010():
  '''
    LEU HB2/HB3 came out swapped (7IT7 fragment). The CH2 reference read the
    CCD's ideal coordinates, and for ARG CB/CG, ILE CG1, LEU CB and MET CB/CG
    those contradict the same entry's model coordinates. The model coordinates
    are what reduce and deposited models follow; all three groups here are -1.
    GLU CB/CG (the CCD agrees with itself) are the control.
  '''
  atoms = place_atoms(pdb_str_010)
  assert hand(atoms, 'LEU', '107', 'CB', 'CA', 'CG', 'HB2') == -1
  assert hand(atoms, 'GLU', '41', 'CB', 'CA', 'CG', 'HB2') == -1
  assert hand(atoms, 'GLU', '41', 'CG', 'CB', 'CD', 'HG2') == -1

def test_011():
  '''
    A propeller came out as the mirror of the CCD every time: the H are still
    superposed when check_propeller_order looks at them, so the order came from
    the riding frame (HB1 got n=2, HB3 n=0). CH3 and NH3 alike; the CCD's two
    coordinate sets agree for 94% of propellers, and here both say -1.
  '''
  atoms = place_atoms(pdb_str_011)
  assert hand(atoms, 'ALA', '42', 'CB', 'HB1', 'HB2', 'CA') == -1
  assert hand(atoms, 'LYS', '87', 'NZ', 'HZ1', 'HZ2', 'CE') == -1
  atoms = place_atoms(pdb_str_010)
  assert hand(atoms, 'LEU', '107', 'CD1', 'HD11', 'HD12', 'CG') == -1
  assert hand(atoms, 'LEU', '107', 'CD2', 'HD21', 'HD22', 'CG') == -1
  # the CH2 of the same residue is unchanged
  assert hand(atoms, 'LEU', '107', 'CB', 'CA', 'CG', 'HB2') == -1

def test_012():
  '''
    AMP (5msd) is a free nucleotide: its O3' is a hydroxyl, but HO3' was
    dropped from every residue that tests as RNA/DNA, where O3' usually
    carries the next phosphate. The H is placed now and removed only where
    the O really is esterified (exclude_H_on_esterified_O).

    The DNA control keeps no HO3' anywhere, including the 3'-terminal DC 3:
    the monomer library entry of a polymer nucleotide has no HO3' to place.
    Giving a 3' end its hydroxyl H is a separate gap, in the library rather
    than here - the RNA entries A/U/G/C have no HO3' either.
  '''
  atoms = place_atoms(pdb_str_012)
  assert h_on(atoms, 'AMP', '1201', "O3'") == ["HO3'"], h_on(
    atoms, 'AMP', '1201', "O3'")
  assert h_on(atoms, 'AMP', '1201', "O2'") == ["HO2'"]
  atoms = place_atoms(pdb_str_012_control)
  for resname, resseq in (('DC', '1'), ('DG', '2'), ('DC', '3')):
    assert h_on(atoms, resname, resseq, "O3'") == [], (resname, resseq)
  # 5GP has the same ligand-style dictionary as AMP, but here its O3' carries
  # the phosphate of C 2 (4TNA), so the H is placed and then removed again
  atoms = place_atoms(pdb_str_012_linked)
  assert h_on(atoms, '5GP', '1', "O3'") == [], h_on(atoms, '5GP', '1', "O3'")
  assert h_on(atoms, '5GP', '1', "O2'") == ["HO2'"]

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

pdb_str_013 = pdb_str_009

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

pdb_str_011 = """
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM     10  N   ALA A  42      12.752 -16.313  13.130  1.00 37.05           N
ATOM     11  CA  ALA A  42      11.861 -16.326  11.950  1.00 41.68           C
ATOM     12  C   ALA A  42      12.606 -15.740  10.747  1.00 42.80           C
ATOM     13  O   ALA A  42      12.597 -16.391   9.684  1.00 36.14           O
ATOM     14  CB  ALA A  42      10.567 -15.603  12.226  1.00 43.61           C
ATOM     61  N   LYS A  87      14.521  55.475  65.158  1.00 16.69           N
ATOM     62  CA  LYS A  87      13.600  54.884  64.201  1.00 17.21           C
ATOM     63  C   LYS A  87      12.210  55.516  64.232  1.00 15.78           C
ATOM     64  O   LYS A  87      11.790  56.116  65.226  1.00 14.71           O
ATOM     65  CB  LYS A  87      13.500  53.368  64.403  1.00 18.88           C
ATOM     66  CG  LYS A  87      13.338  52.939  65.839  1.00 24.38           C
ATOM     67  CD  LYS A  87      11.961  52.424  66.104  1.00 30.22           C
ATOM     68  CE  LYS A  87      11.700  51.121  65.394  1.00 26.16           C
ATOM     69  NZ  LYS A  87      10.264  50.804  65.531  1.00 31.31           N
END
"""

pdb_str_012 = """
CRYST1   52.650   54.499   66.365  93.37  97.68 102.44 P 1
HETATM    1  C1' AMP A1201       1.751   9.985  45.277  1.00 24.11           C
HETATM    2  C2  AMP A1201       3.845   8.243  48.752  1.00 26.36           C
HETATM    3  C2' AMP A1201       2.736  11.117  44.952  1.00 24.56           C
HETATM    4  C3' AMP A1201       2.584  11.238  43.440  1.00 24.53           C
HETATM    5  C4  AMP A1201       3.059   8.342  46.558  1.00 23.44           C
HETATM    6  C4' AMP A1201       1.117  10.895  43.217  1.00 26.37           C
HETATM    7  C5  AMP A1201       3.693   7.052  46.280  1.00 24.64           C
HETATM    8  C5' AMP A1201       0.762  10.378  41.845  1.00 27.37           C
HETATM    9  C6  AMP A1201       4.423   6.413  47.392  1.00 24.61           C
HETATM   10  C8  AMP A1201       2.694   7.786  44.489  1.00 23.99           C
HETATM   11  N1  AMP A1201       4.457   7.065  48.572  1.00 22.89           N
HETATM   12  N3  AMP A1201       3.166   8.890  47.789  1.00 25.59           N
HETATM   13  N6  AMP A1201       5.046   5.216  47.246  1.00 26.07           N
HETATM   14  N7  AMP A1201       3.442   6.772  44.990  1.00 23.32           N
HETATM   15  N9  AMP A1201       2.467   8.713  45.433  1.00 24.84           N
HETATM   16  O1P AMP A1201      -0.476  10.260  39.151  1.00 31.62           O
HETATM   17  O2' AMP A1201       2.299  12.295  45.610  1.00 27.03           O
HETATM   18  O2P AMP A1201      -1.283   8.318  40.533  1.00 33.06           O
HETATM   19  O3' AMP A1201       2.930  12.509  42.922  1.00 26.18           O
HETATM   20  O3P AMP A1201      -2.715  10.330  40.359  1.00 28.38           O
HETATM   21  O4' AMP A1201       0.851   9.868  44.194  1.00 25.17           O
HETATM   22  O5' AMP A1201      -0.643  10.350  41.667  1.00 29.44           O
HETATM   23  P   AMP A1201      -1.307   9.798  40.318  1.00 33.89           P
END
"""

pdb_str_012_control = """
CRYST1   17.880   31.420   43.900  90.00  90.00  90.00 P 21 21 21    8
ATOM      1  O5'  DC A   1      19.545  18.136  17.917  1.00  3.07           O
ATOM      2  C5'  DC A   1      19.769  17.119  18.884  1.00  2.46           C
ATOM      3  C4'  DC A   1      18.610  16.148  19.001  1.00  2.19           C
ATOM      4  O4'  DC A   1      17.462  16.852  19.514  1.00  2.62           O
ATOM      5  C3'  DC A   1      18.161  15.506  17.674  1.00  2.26           C
ATOM      6  O3'  DC A   1      17.782  14.139  17.875  1.00  2.47           O
ATOM      7  C2'  DC A   1      16.906  16.282  17.315  1.00  2.41           C
ATOM      8  C1'  DC A   1      16.340  16.624  18.692  1.00  2.34           C
ATOM      9  N1   DC A   1      15.516  17.837  18.704  1.00  2.30           N
ATOM     10  C2   DC A   1      14.145  17.720  18.492  1.00  2.34           C
ATOM     11  O2   DC A   1      13.658  16.581  18.329  1.00  3.07           O
ATOM     12  N3   DC A   1      13.385  18.831  18.454  1.00  2.35           N
ATOM     13  C4   DC A   1      13.943  20.039  18.611  1.00  2.37           C
ATOM     14  N4   DC A   1      13.161  21.108  18.580  1.00  2.79           N
ATOM     15  C5   DC A   1      15.357  20.189  18.812  1.00  2.71           C
ATOM     16  C6   DC A   1      16.103  19.061  18.841  1.00  2.68           C
ATOM     17  P    DG A   2      18.825  12.942  17.684  1.00  2.51           P
ATOM     18  OP1  DG A   2      19.788  13.206  16.573  1.00  3.24           O
ATOM     19  OP2  DG A   2      17.976  11.710  17.621  1.00  3.25           O
ATOM     20  O5'  DG A   2      19.719  12.937  19.002  1.00  2.66           O
ATOM     21  C5'  DG A   2      19.103  12.733  20.284  1.00  2.85           C
ATOM     22  C4'  DG A   2      20.140  13.045  21.335  1.00  2.73           C
ATOM     23  O4'  DG A   2      20.546  14.399  21.207  1.00  2.76           O
ATOM     24  C3'  DG A   2      19.598  12.910  22.753  1.00  2.85           C
ATOM     25  O3'  DG A   2      19.812  11.563  23.230  1.00  3.51           O
ATOM     26  C2'  DG A   2      20.430  13.919  23.526  1.00  3.30           C
ATOM     27  C1'  DG A   2      20.834  14.964  22.481  1.00  2.85           C
ATOM     28  N9   DG A   2      20.140  16.222  22.572  1.00  2.76           N
ATOM     29  C8   DG A   2      20.744  17.451  22.654  1.00  3.56           C
ATOM     30  N7   DG A   2      19.903  18.441  22.626  1.00  3.75           N
ATOM     31  C5   DG A   2      18.658  17.840  22.510  1.00  2.68           C
ATOM     32  C6   DG A   2      17.359  18.399  22.397  1.00  2.62           C
ATOM     33  O6   DG A   2      17.067  19.612  22.382  1.00  3.52           O
ATOM     34  N1   DG A   2      16.371  17.435  22.313  1.00  2.19           N
ATOM     35  C2   DG A   2      16.610  16.084  22.263  1.00  1.99           C
ATOM     36  N2   DG A   2      15.540  15.294  22.108  1.00  2.48           N
ATOM     37  N3   DG A   2      17.814  15.534  22.356  1.00  2.10           N
ATOM     38  C4   DG A   2      18.782  16.466  22.477  1.00  2.26           C
ATOM     39  P    DC A   3      18.598  10.645  23.664  1.00  3.54           P
ATOM     40  OP1  DC A   3      19.187   9.314  23.999  1.00  5.42           O
ATOM     41  OP2  DC A   3      17.526  10.682  22.656  1.00  5.19           O
ATOM     42  O5'  DC A   3      18.133  11.412  24.998  1.00  3.68           O
ATOM     43  C5'  DC A   3      17.562  10.688  26.104  1.00  3.12           C
ATOM     44  C4'  DC A   3      16.175  11.192  26.422  1.00  2.67           C
ATOM     45  O4'  DC A   3      16.263  12.557  26.885  1.00  2.93           O
ATOM     46  C3'  DC A   3      15.179  11.206  25.269  1.00  2.83           C
ATOM     47  O3'  DC A   3      13.854  10.971  25.771  1.00  3.13           O
ATOM     48  C2'  DC A   3      15.208  12.650  24.803  1.00  3.07           C
ATOM     49  C1'  DC A   3      15.379  13.367  26.137  1.00  2.64           C
ATOM     50  N1   DC A   3      15.979  14.696  26.016  1.00  2.59           N
ATOM     51  C2   DC A   3      15.126  15.784  25.823  1.00  2.60           C
ATOM     52  O2   DC A   3      13.894  15.575  25.734  1.00  3.08           O
ATOM     53  N3   DC A   3      15.658  17.016  25.725  1.00  2.61           N
ATOM     54  C4   DC A   3      16.987  17.184  25.770  1.00  2.65           C
ATOM     55  N4   DC A   3      17.464  18.422  25.674  1.00  3.18           N
ATOM     56  C5   DC A   3      17.881  16.072  25.908  1.00  3.00           C
ATOM     57  C6   DC A   3      17.330  14.855  26.026  1.00  2.98           C
END
"""

pdb_str_012_linked = """
CRYST1   56.300   33.400   63.000  90.00  90.25  90.00 P 1 21 1      2
ATOM      1  OP3 5GP A   1      23.215   5.145  51.161  1.00  0.00           O
ATOM      2  P   5GP A   1      24.650   4.594  51.620  1.00  0.00           P
ATOM      3  OP1 5GP A   1      24.810   3.219  51.082  1.00  0.00           O
ATOM      4  OP2 5GP A   1      25.707   5.598  51.344  1.00  0.00           O
ATOM      5  O5' 5GP A   1      24.476   4.546  53.204  1.00  0.00           O
ATOM      6  C5' 5GP A   1      25.156   5.525  53.978  1.00  0.00           C
ATOM      7  C4' 5GP A   1      25.808   4.801  55.152  1.00  0.00           C
ATOM      8  O4' 5GP A   1      24.983   3.707  55.565  1.00  0.00           O
ATOM      9  C3' 5GP A   1      27.185   4.241  54.785  1.00  0.00           C
ATOM     10  O3' 5GP A   1      28.130   5.240  55.163  1.00  0.00           O
ATOM     11  C2' 5GP A   1      27.282   3.124  55.808  1.00  0.00           C
ATOM     12  O2' 5GP A   1      27.691   3.684  57.045  1.00  0.00           O
ATOM     13  C1' 5GP A   1      25.833   2.612  55.903  1.00  0.00           C
ATOM     14  N9  5GP A   1      25.681   1.623  54.832  1.00  0.00           N
ATOM     15  C8  5GP A   1      24.865   1.712  53.755  1.00  0.00           C
ATOM     16  N7  5GP A   1      25.012   0.580  52.989  1.00  0.00           N
ATOM     17  C5  5GP A   1      25.874  -0.248  53.697  1.00  0.00           C
ATOM     18  C6  5GP A   1      26.336  -1.541  53.483  1.00  0.00           C
ATOM     19  O6  5GP A   1      25.983  -2.159  52.466  1.00  0.00           O
ATOM     20  N1  5GP A   1      27.127  -2.122  54.405  1.00  0.00           N
ATOM     21  C2  5GP A   1      27.481  -1.468  55.532  1.00  0.00           C
ATOM     22  N2  5GP A   1      28.055  -2.114  56.534  1.00  0.00           N
ATOM     23  N3  5GP A   1      27.078  -0.243  55.747  1.00  0.00           N
ATOM     24  C4  5GP A   1      26.260   0.388  54.839  1.00  0.00           C
ATOM     25  P     C A   2      29.640   5.301  54.591  1.00  0.00           P
ATOM     26  OP1   C A   2      30.307   6.511  55.152  1.00  0.00           O
ATOM     27  OP2   C A   2      29.605   5.142  53.117  1.00  0.00           O
ATOM     28  O5'   C A   2      30.332   4.036  55.292  1.00  0.00           O
ATOM     29  C5'   C A   2      30.763   4.137  56.639  1.00  0.00           C
ATOM     30  C4'   C A   2      31.424   2.812  56.982  1.00  0.00           C
ATOM     31  O4'   C A   2      30.483   1.753  56.828  1.00  0.00           O
ATOM     32  C3'   C A   2      32.562   2.541  56.002  1.00  0.00           C
ATOM     33  O3'   C A   2      33.772   3.124  56.479  1.00  0.00           O
ATOM     34  C2'   C A   2      32.637   1.038  56.077  1.00  0.00           C
ATOM     35  O2'   C A   2      33.390   0.643  57.211  1.00  0.00           O
ATOM     36  C1'   C A   2      31.162   0.636  56.248  1.00  0.00           C
ATOM     37  N1    C A   2      30.570   0.343  54.911  1.00  0.00           N
ATOM     38  C2    C A   2      30.799  -0.845  54.338  1.00  0.00           C
ATOM     39  O2    C A   2      31.777  -1.511  54.728  1.00  0.00           O
ATOM     40  N3    C A   2      30.243  -1.143  53.141  1.00  0.00           N
ATOM     41  C4    C A   2      29.496  -0.265  52.481  1.00  0.00           C
ATOM     42  N4    C A   2      28.782  -0.652  51.425  1.00  0.00           N
ATOM     43  C5    C A   2      29.284   0.988  53.032  1.00  0.00           C
ATOM     44  C6    C A   2      29.844   1.264  54.273  1.00  0.00           C
END
"""

if __name__ == '__main__':
  run()
  print("OK")
