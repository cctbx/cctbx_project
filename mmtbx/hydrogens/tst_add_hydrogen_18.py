"""place_hydrogens leaves the received model unchanged: hierarchy, restraint
object list (same object, same entries), restraints manager and crystal
symmetry, for keep_existing_H True and False, with a cell, without symmetry
and with an element X atom. The list is shared with the model's copies
(select(), deep_copy()), so a CCD-only ligand's auto_ dictionary must not be
appended to it."""
from __future__ import absolute_import, division, print_function
import iotbx.cif
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

cryst1 = "CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1\n"

# 1yjp residues 1-3, ligand 04 (user CIF), BGM (CCD only, no geostd entry)
pdb_str = """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
HETATM   22  C1   04 B   1      15.000  10.000  15.000  1.00 20.00           C
HETATM   23  C2   04 B   1      16.520  10.000  15.000  1.00 20.00           C
HETATM   24  O1   04 B   1      17.220  11.210  15.000  1.00 20.00           O
HETATM   25  P   BGM C   1      -9.873  19.915  28.793  1.00 20.00           P
HETATM   26  OP1 BGM C   1     -10.967  19.590  29.736  1.00 20.00           O
HETATM   27  OP2 BGM C   1      -9.665  21.510  28.742  1.00 20.00           O
HETATM   28  O5' BGM C   1      -8.512  19.211  29.289  1.00 20.00           O
HETATM   29  C5' BGM C   1      -7.969  19.403  30.597  1.00 20.00           C
HETATM   30  C4' BGM C   1      -6.677  18.594  30.735  1.00 20.00           C
HETATM   31  O4' BGM C   1      -5.653  19.118  29.863  1.00 20.00           O
HETATM   32  C1' BGM C   1      -4.389  18.682  30.389  1.00 20.00           C
HETATM   33  N9  BGM C   1      -3.370  19.708  30.154  1.00 20.00           N
HETATM   34  C8  BGM C   1      -3.522  21.053  30.334  1.00 20.00           C
HETATM   35  N7  BGM C   1      -2.414  21.666  30.030  1.00 20.00           N
HETATM   36  C5  BGM C   1      -1.490  20.754  29.638  1.00 20.00           C
HETATM   37  C4  BGM C   1      -2.090  19.493  29.718  1.00 20.00           C
HETATM   38  N3  BGM C   1      -1.388  18.403  29.389  1.00 20.00           N
HETATM   39  C2  BGM C   1      -0.138  18.494  28.992  1.00 20.00           C
HETATM   40  N2  BGM C   1       0.542  17.348  28.664  1.00 20.00           N
HETATM   41  N1  BGM C   1       0.501  19.694  28.896  1.00 20.00           N
HETATM   42  C6  BGM C   1      -0.147  20.839  29.205  1.00 20.00           C
HETATM   43  O6  BGM C   1       0.417  21.917  29.124  1.00 20.00           O
HETATM   44  C2' BGM C   1      -4.574  18.465  31.904  1.00 20.00           C
HETATM   45  C3' BGM C   1      -6.071  18.772  32.151  1.00 20.00           C
HETATM   46  O3' BGM C   1      -6.639  17.837  33.070  1.00 20.00           O
HETATM   47  OP3 BGM C   1     -10.252  19.375  27.325  1.00 20.00           O
HETATM   48  BR  BGM C   1      -5.103  21.898  30.935  1.00 20.00          BR
"""

# ethanol, as a user dictionary under a code that is not in the CCD
cif_04 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 04  04  'test ligand'  ligand  9  3  .

data_comp_04
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
 04  C1   C  CH3   0
 04  C2   C  CH2   0
 04  O1   O  OH1   0
 04  H11  H  HCH3  0
 04  H12  H  HCH3  0
 04  H13  H  HCH3  0
 04  H21  H  HCH2  0
 04  H22  H  HCH2  0
 04  HO1  H  HOH1  0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 04  C1  C2   single  1.520  0.020
 04  C2  O1   single  1.430  0.020
 04  C1  H11  single  0.970  0.020
 04  C1  H12  single  0.970  0.020
 04  C1  H13  single  0.970  0.020
 04  C2  H21  single  0.970  0.020
 04  C2  H22  single  0.970  0.020
 04  O1  HO1  single  0.850  0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 04  C1   C2  O1   109.5  3.0
 04  C2   C1  H11  109.5  3.0
 04  C2   C1  H12  109.5  3.0
 04  C2   C1  H13  109.5  3.0
 04  H11  C1  H12  109.5  3.0
 04  H11  C1  H13  109.5  3.0
 04  H12  C1  H13  109.5  3.0
 04  C1   C2  H21  109.5  3.0
 04  C1   C2  H22  109.5  3.0
 04  O1   C2  H21  109.5  3.0
 04  O1   C2  H22  109.5  3.0
 04  H21  C2  H22  109.5  3.0
 04  C2   O1  HO1  109.5  3.0
"""

element_x = (
  "HETATM   49  X   UNX D   1      10.000  10.000  10.000  1.00 20.00"
  "           X\n")

def caller_model(text, processed):
  ro = [("04.cif", iotbx.cif.reader(input_string=cif_04).model())]
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=text),
    restraint_objects=ro, stop_for_unknowns=False, log=null_out())
  if processed:
    model.process(make_restraints=True)
  return model

def snapshot(model):
  ro = model.get_restraint_objects()
  return dict(
    atoms=[(a.id_str(), a.xyz) for a in model.get_hierarchy().atoms()],
    hierarchy=model.get_hierarchy(),
    ro=ro, ro_entries=list(ro),
    rm=model.get_restraints_manager(),
    cs=model.crystal_symmetry())

def check_unchanged(model, before, label):
  after = snapshot(model)
  assert after["atoms"] == before["atoms"], label
  assert after["hierarchy"] is before["hierarchy"], label
  assert after["ro"] is before["ro"], label
  assert after["ro_entries"] == before["ro_entries"], label
  assert [n for n, o in after["ro"]] == ["04.cif"], (label, after["ro"])
  assert after["rm"] is before["rm"], label
  if before["cs"] is None:
    assert after["cs"] is None, label
  else:
    assert after["cs"].is_identical_symmetry(before["cs"]), label

def h_names(model, resname):
  return sorted(a.name.strip() for a in model.get_hierarchy().atoms()
    if a.parent().resname.strip() == resname and a.element.strip() == "H")

def run_case(text, processed, label, with_copy=True):
  for keep in (False, True):
    lab = "%s, keep_existing_H=%s" % (label, keep)
    model = caller_model(text, processed)
    # deep_copy() stops on element X (no scattering type)
    copy = model.deep_copy() if with_copy else None
    if copy is not None:
      assert copy.get_restraint_objects() is model.get_restraint_objects()
    before = snapshot(model)
    obj = reduce_hydrogen.place_hydrogens(model=model, keep_existing_H=keep)
    obj.run()
    placed = obj.get_model()
    check_unchanged(model, before, lab)
    if copy is not None:
      assert [n for n, o in copy.get_restraint_objects()] == ["04.cif"], lab
    # the placed model still gets H on BGM (CCD only) and on 04 (user CIF)
    assert len(h_names(placed, "BGM")) == 13, (lab, h_names(placed, "BGM"))
    assert h_names(placed, "04") == \
      ["H11", "H12", "H13", "H21", "H22", "HO1"], (lab, h_names(placed, "04"))
    assert placed.get_hierarchy() is not model.get_hierarchy(), lab
    assert not [n for n, o in placed.get_restraint_objects()
                if n.startswith("auto_")], lab

def test_with_cell():
  run_case(cryst1 + pdb_str, processed=True, label="cell, processed")
  run_case(cryst1 + pdb_str, processed=False, label="cell, unprocessed")

def test_without_symmetry():
  run_case(pdb_str, processed=False, label="no symmetry, unprocessed")
  run_case(pdb_str, processed=True, label="no symmetry, processed")

def test_element_x():
  run_case(cryst1 + pdb_str + element_x, processed=False,
    label="element X, cell", with_copy=False)
  run_case(pdb_str + element_x, processed=False,
    label="element X, no symmetry", with_copy=False)

def test_add_auto_restraint_object():
  """the model gets a new list, its copies keep the old one"""
  model = caller_model(cryst1 + pdb_str, processed=False)
  copy = model.deep_copy()
  ro = model.get_restraint_objects()
  names = []
  reduce_hydrogen.add_auto_restraint_object(model, "BGM", "cif", names)
  assert names == ["auto_BGM"]
  assert [n for n, o in ro] == ["04.cif"]
  assert copy.get_restraint_objects() is ro
  assert [n for n, o in model.get_restraint_objects()] == \
    ["04.cif", "auto_BGM"]
  reduce_hydrogen.remove_auto_restraint_objects(model, names)
  assert [n for n, o in model.get_restraint_objects()] == ["04.cif"]

if (__name__ == "__main__"):
  test_with_cell()
  test_without_symmetry()
  test_element_x()
  test_add_auto_restraint_object()
  print("OK")
