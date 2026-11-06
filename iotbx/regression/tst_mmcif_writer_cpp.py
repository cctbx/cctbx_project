"""
The mmCIF writer, iotbx.pdb.hierarchy.root.as_cif_block, against its Python
body, root._as_cif_block_python, kept as the reference while the atom loop
moves to C++: the two give the same block for every fixture below. Every
use of the reference is in a block marked PARITY.
"""
import os
import libtbx.load_env
from cctbx.array_family import flex
import iotbx.pdb
import iotbx.pdb.hierarchy
import iotbx.pdb.tst_mmcif_hierarchy as tst_mmcif_hierarchy


# Two models (the serials restart), ANISOU in both, altlocs with an
# insertion code, a BREAK, two charged atoms, one a HETATM ion, and a water
# chain without segid.
multi_model_pdb = """\
CRYST1   20.000   30.000   40.000  90.00  90.00  90.00 P 1
MODEL        1
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00      SEGA N
ANISOU    1  N   ALA A   1      100    200    300     10     20     30  SEGA N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00 11.00      SEGA C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00 12.00      SEGA C
ATOM      4  N  AGLY A   2A      4.000   2.000   3.000  0.50 13.00      SEGA N
ATOM      5  N  BGLY A   2A      4.100   2.100   3.100  0.50 13.00      SEGA N
BREAK
ATOM      6  N   SER A   5       9.000   2.000   3.000  1.00 14.00      SEGA N1+
HETATM    7 ZN    ZN A 101      10.000   2.000   3.000  1.00 15.00      SEGAZN2+
TER
HETATM    8  O   HOH B   1      11.000   2.000   3.000  1.00 16.00           O
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1       1.500   2.000   3.000  1.00 10.00      SEGA N
ANISOU    1  N   ALA A   1      110    210    310     11     21     31  SEGA N
ATOM      2  CA  ALA A   1       2.500   2.000   3.000  1.00 11.00      SEGA C
ENDMDL
END
"""


# 27 zinc ions in one chain: each a ligand with its own label_asym_id, more
# ids than letters, so that the first-appearance order of _struct_asym
# differs from the sorted one (AA sorts before B)
many_ligands_pdb = "".join(
  "HETATM%5d ZN    ZN A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          ZN\n"
  % (i + 1, i + 1, 3.0 * i, 0, 0, 1, 10) for i in range(27))


# The model files of iotbx/regression/data, named so that a file added there
# for another test cannot break this one or its time budget. The .pdb twins
# of the three non_zero_origin .cif files are left out: their blocks are
# the same.
data_files = (
  "1yjp.pdb", "2ERL.pdb", "big_cube_model.cif", "big_cube_model.pdb",
  "non_zero_origin_model.cif", "non_zero_origin_model_split.cif",
  "non_zero_origin_ncs_model.cif", "parity_3orl.cif", "parity_ca_trace.cif",
  "parity_icode_altloc.cif", "parity_ligand_only.cif",
  "parity_protein_ligand.cif", "parity_two_models.cif", "pdb1zff.ent")


def fixtures():
  """The inputs of the parity checks.

  Returns
  -------
  list of tuple
    (name, iotbx.pdb.hierarchy.root, cctbx.crystal.symmetry or None) for
    six inputs of iotbx/pdb/tst_mmcif_hierarchy.py, multi_model_pdb,
    many_ligands_pdb and the data_files.
  """
  t = tst_mmcif_hierarchy
  # the inputs of tst_mmcif_hierarchy.py with a feature the other fixtures
  # lack: no break records (input_1ab1), anisotrop and symmetry from mmCIF
  # (input_4edr), fp and fdp with charges (input_anom), a four-character
  # residue number (input_charges), D-amino acids with waters and altlocs
  # (input_4gln), no atoms (atom_site_header)
  texts = [
    ("input_1ab1", t.input_1ab1),
    ("input_4edr", t.input_4edr),
    ("input_anom", t.input_anom),
    ("input_charges", t.input_charges),
    ("input_4gln", t.pdb_atom_site_loop_header + t.input_4gln),
    ("atom_site_header", t.atom_site_header),
    ("multi_model_pdb", multi_model_pdb),
    ("many_ligands_pdb", many_ligands_pdb),
  ]
  result = []
  for name, text in texts:
    inp = iotbx.pdb.input(source_info=None, lines=flex.split_lines(text))
    result.append((name, inp.construct_hierarchy(), inp.crystal_symmetry()))
  data_dir = libtbx.env.under_dist("iotbx", "regression/data")
  for name in data_files:
    inp = iotbx.pdb.input(file_name=os.path.join(data_dir, name))
    result.append((name, inp.construct_hierarchy(), inp.crystal_symmetry()))
  return result


def fixture(fixture_list, name):
  """The (hierarchy, crystal symmetry) of the fixture called name."""
  for fixture_name, h, cs in fixture_list:
    if fixture_name == name:
      return h, cs
  raise KeyError(name)


def compare(name, h, cs, **kw):
  """The same block from as_cif_block and the Python reference."""
  # PARITY-BEGIN
  from libtbx.test_utils import show_diff
  got = str(h.as_cif_block(cs, **kw))
  expected = str(h._as_cif_block_python(cs, **kw))
  assert not show_diff(got, expected), (name, kw)
  # PARITY-END


def build_root(model_ids, serial=None):
  """One chain A with one ALA residue of two atoms in each model; default
  atom fields apart from the name and, when given, the serial."""
  root = iotbx.pdb.hierarchy.root()
  for model_id in model_ids:
    model = iotbx.pdb.hierarchy.model(id=model_id)
    root.append_model(model)
    chain = iotbx.pdb.hierarchy.chain(id="A")
    model.append_chain(chain)
    rg = iotbx.pdb.hierarchy.residue_group(resseq="   1")
    chain.append_residue_group(rg)
    ag = iotbx.pdb.hierarchy.atom_group(resname="ALA")
    rg.append_atom_group(ag)
    for atom_name in (" N  ", " CA "):
      atom = iotbx.pdb.hierarchy.atom()
      atom.set_name(atom_name)
      if serial is not None:
        atom.serial = serial
      ag.append_atom(atom)
  return root


def exercise_reference_exists(fixture_list):
  for name, h, cs in fixture_list:
    compare(name, h, cs)


def exercise_atom_site_labels(fixture_list):
  h, cs = fixture(fixture_list, "multi_model_pdb")
  # model 1: chain A (ALA 1, GLY 2A in two conformers, SER 5 after the
  # BREAK, ZN 101), chain B (HOH 1); model 2: chain A (ALA 1). The segid
  # SEGA is longer than the chain id, so it is the auth_asym_id of A unless
  # it is written as auth_segid.
  for segid_as_auth_segid, chain_auth_asym_id in (
      (False, ["SEGA", "B", "SEGA"]), (True, ["A", "B", "A"])):
    labels = h._atom_site_labels(segid_as_auth_segid=segid_as_auth_segid)
    assert list(labels.chain_auth_asym_id) == chain_auth_asym_id
    assert list(labels.rg_label_asym_id) == ["A", "A", "A", "B", "C", "E"]
    assert list(labels.rg_auth_seq_id) == ["1", "2", "5", "101", "1", "1"]
    assert list(labels.ag_label_seq_id) == [
      "1", "2", "2", "3", ".", ".", "1"]
  # a hybrid-36 residue number is decoded
  root = build_root([""], serial="    1")
  root.only_residue_group().resseq = "A000"
  assert list(root._atom_site_labels().rg_auth_seq_id) == ["10000"]
  # the lookups are rebuilt on every call: a chain added after a call gets
  # its label_asym_id
  h = h.deep_copy()
  h._atom_site_labels()
  atom = iotbx.pdb.hierarchy.atom()
  atom.set_name(" O  ")
  ag = iotbx.pdb.hierarchy.atom_group(resname="HOH")
  ag.append_atom(atom)
  rg = iotbx.pdb.hierarchy.residue_group(resseq="   7")
  rg.append_atom_group(ag)
  chain = iotbx.pdb.hierarchy.chain(id="C")
  chain.append_residue_group(rg)
  h.models()[1].append_chain(chain)
  labels = h._atom_site_labels()
  assert list(labels.chain_auth_asym_id) == ["SEGA", "B", "SEGA", "C"]
  assert list(labels.rg_label_asym_id) == ["A", "A", "A", "B", "C", "E", "F"]
  # model 2 now has a chain model 1 lacks: the label arrays run across
  # models
  compare("added chain", h, None)


def run():
  # all_label_asym_ids() builds its list of 475254 ids, A to ZZZZ, on every
  # call, 0.2 s, and every as_cif_block call makes one call; the checks here
  # make dozens, so they share one list, up to ZZZ, more than any of them
  # needs.
  all_label_asym_ids = iotbx.pdb.hierarchy.all_label_asym_ids
  label_asym_ids = all_label_asym_ids(maximum_length=3)
  iotbx.pdb.hierarchy.all_label_asym_ids = lambda: label_asym_ids
  try:
    fixture_list = fixtures()
    exercise_reference_exists(fixture_list)
    exercise_atom_site_labels(fixture_list)
  finally:
    iotbx.pdb.hierarchy.all_label_asym_ids = all_label_asym_ids


if __name__ == "__main__":
  run()
  print("OK")
