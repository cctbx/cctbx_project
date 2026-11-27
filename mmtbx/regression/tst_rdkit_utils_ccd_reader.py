"""
rdkit_utils.read_chemical_component_filename(): which coordinates the atoms
of a CCD entry get. The ideal coordinates if every atom has all three; the
model coordinates if no atom has an ideal coordinate and every atom has all
three model coordinates; else a Sorry that asks for restraints. Never one
atom's coordinates on another atom (MDE: two H without ideal coordinates),
and no model coordinates in place of partly missing ideal ones.
"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
from mmtbx.chemical_components import get_cif_filename
from mmtbx.ligands import rdkit_utils

# a water under the name ZZW; model and ideal coordinates differ
head = """\
data_ZZW
#
_chem_comp.id                                   ZZW
_chem_comp.name                                 WATER
_chem_comp.type                                 NON-POLYMER
_chem_comp.three_letter_code                    ZZW
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
"""
bonds = """\
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
ZZW  O  H1  SING
ZZW  O  H2  SING
"""
model = {"O": (-23.107, 18.401, -21.626), "H1": (-22.157, 18.401, -21.626),
         "H2": (-23.424, 18.401, -20.730)}
ideal = {"O": (-0.064, 0.0, 0.0), "H1": (0.512, 0.0, -0.776),
         "H2": (0.512, 0.0, 0.776)}

def entry(no_model=(), no_ideal=(), mark="?", ideal_marks=None):
  """no_model, no_ideal: atom names without that set (mark in its three
  fields); ideal_marks: {atom name: (x, y, z) with None for the value}"""
  rows = []
  for name in ("O", "H1", "H2"):
    m = [mark] * 3 if name in no_model else ["%.3f" % x for x in model[name]]
    i = [mark] * 3 if name in no_ideal else ["%.3f" % x for x in ideal[name]]
    for n, v in enumerate((ideal_marks or {}).get(name, ())):
      if v is not None: i[n] = v
    rows.append("ZZW  %s  %s  0  %s  %s" % (name, name[0], "  ".join(m),
      "  ".join(i)))
  return head + "\n".join(rows) + "\n" + bonds

def sites(text, prefix):
  """{atom name: xyz} of the molecule read from the entry text"""
  fn = prefix + ".cif"
  with open(fn, "w") as f:
    f.write(text)
  return sites_of(rdkit_utils.read_chemical_component_filename(fn))

def sites_of(mol):
  conformer = mol.GetConformer()
  result = {}
  for atom in mol.GetAtoms():
    p = conformer.GetAtomPosition(atom.GetIdx())
    result[atom.GetProp("atom_id")] = (p.x, p.y, p.z)
  return result

def check(found, expected):
  assert sorted(found) == sorted(expected), sorted(found)
  for name in expected:
    assert approx_equal(found[name], expected[name], eps=1.e-3, out=None), (
      name, found[name], expected[name])

def sorry(text, prefix):
  try:
    sites(text, prefix)
  except Sorry as e:
    assert "Please supply" in str(e), str(e)
    return str(e)
  raise AssertionError("no Sorry: %s" % prefix)

def exercise_all_ideal():
  check(sites(entry(), "tst_ccd_reader_1"), ideal)

def exercise_no_ideal():
  # no ideal coordinates at all: the complete model coordinates, as before
  for i, mark in enumerate("?."):
    check(sites(entry(no_ideal=("O", "H1", "H2"), mark=mark),
      "tst_ccd_reader_2_%d" % i), model)

def exercise_some_ideal_missing():
  # was: H2 dropped from the list, IndexError at the last atom; with the gap
  # at H1, H2's ideal coordinates went to H1 first. Now a Sorry that names
  # the atoms; the complete model coordinates are not used instead.
  text = sorry(entry(no_ideal=("H2",)), "tst_ccd_reader_3")
  assert "for ZZW failed: 1 of its 3 atoms has no ideal coordinates (H2)." \
    in text, text
  assert "\n" not in text.strip(), text
  text = sorry(entry(no_ideal=("H1", "H2"), mark="."), "tst_ccd_reader_4")
  assert "2 of its 3 atoms have no ideal coordinates (H1 H2)." in text, text

def exercise_one_component_missing():
  # one missing component is a missing row (was: H1 at 0, 0, 0 for x = ?;
  # ValueError for y = ? or for '.')
  for i, marks in enumerate((("?", None, None), (None, "?", None),
                             (None, None, "."), (".", None, "?"))):
    text = sorry(entry(ideal_marks={"H1": marks}),
      "tst_ccd_reader_6_%d" % i)
    assert "1 of its 3 atoms has incomplete ideal coordinates (H1)." in \
      text, text

def exercise_one_component_missing_everywhere():
  # ideal x unknown on every atom, y and z given: the set and the atoms are
  # named (was: the text of an entry without any coordinates)
  marks = dict((name, ("?", None, None)) for name in ("O", "H1", "H2"))
  text = sorry(entry(ideal_marks=marks), "tst_ccd_reader_7")
  assert ("for ZZW failed: 3 of its 3 atoms have incomplete ideal "
    "coordinates (O H1 H2). Please supply restraints.") in text, text

def exercise_restraints_file():
  # read_restraints_file(): the plain x, y, z of a restraints file
  text = entry(no_ideal=("H2",))
  for x in "xyz":
    text = text.replace("_chem_comp_atom.pdbx_model_Cartn_%s_ideal" % x,
      "_chem_comp_atom.%s" % x)
  fn = "tst_ccd_reader_8.cif"
  with open(fn, "w") as f:
    f.write(text)
  try:
    rdkit_utils.read_restraints_file(fn)
  except Sorry as e:
    assert "1 of its 3 atoms has no x, y, z coordinates (H2)." in str(e), \
      str(e)
    assert "ideal" not in str(e), str(e)
  else:
    raise AssertionError("no Sorry")
  with open(fn, "w") as f:
    f.write(text.replace("?  ?  ?", "0.512  0.000  0.776"))
  check(sites_of(rdkit_utils.read_restraints_file(fn)), ideal)

def exercise_model_incomplete():
  # no ideal coordinates, and the model coordinates have a gap (was:
  # IndexError); none of either
  text = sorry(entry(no_ideal=("O", "H1", "H2"), no_model=("H2",)),
    "tst_ccd_reader_5_0")
  assert "1 of its 3 atoms has no model coordinates (H2)" in text, text
  sorry(entry(no_ideal=("O", "H1", "H2"), no_model=("O", "H1", "H2")),
    "tst_ccd_reader_5_1")
  # ideal coordinates with a gap stop whatever the model coordinates are
  text = sorry(entry(no_ideal=("H1",), no_model=("H2",)),
    "tst_ccd_reader_5_2")
  assert "no ideal coordinates (H1)" in text, text

def exercise_mde():
  # CCD MDE: HO2 and HO5 have no ideal coordinates
  fn = get_cif_filename("MDE")
  if not os.path.isfile(fn):
    print("Skipping exercise_mde: no CCD"); return
  try:
    rdkit_utils.read_chemical_component_filename(fn)
  except Sorry as e:
    assert ("for MDE failed: 2 of its 118 atoms have no ideal coordinates "
      "(HO2 HO5). Please supply restraints.") in str(e), str(e)
  else:
    raise AssertionError("no Sorry")
  assert rdkit_utils.get_molecule_from_resname("MDE") is None

if __name__ == "__main__":
  exercise_all_ideal()
  exercise_no_ideal()
  exercise_some_ideal_missing()
  exercise_one_component_missing()
  exercise_one_component_missing_everywhere()
  exercise_restraints_file()
  exercise_model_incomplete()
  exercise_mde()
  print("OK")
