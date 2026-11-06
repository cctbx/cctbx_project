"""
The mmCIF writer, iotbx.pdb.hierarchy.root.as_cif_block: the block holds
the C++ columns (iotbx.pdb.atom_site_columns_from_hierarchy) under their
tags (check_block), the columns hold the values the exercises derive by
hand, and both equal those of the Python body as_cif_block replaced,
root._as_cif_block_python, for every fixture below. Every use of that
reference is in a block marked PARITY. The two differences from it are
pinned: label_alt_id of an atom whose parent() is not the atom group that
holds it (exercise_atoms_with_other_parent), and the exception of an input
with two defects (exercise_errors). Retiring the reference deletes the
PARITY blocks; compare and compare_columns then call check_block and
cpp_columns directly in place of themselves, and the comments that name
the reference are reworded.
"""
import locale
import os
import libtbx.load_env
from libtbx.test_utils import Exception_expected
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


# The atom_site_output members and the tags of their columns, in the order
# of the loop.
atom_site_tags = (
  ("group_pdb", "_atom_site.group_PDB"),
  ("id", "_atom_site.id"),
  ("label_atom_id", "_atom_site.label_atom_id"),
  ("label_alt_id", "_atom_site.label_alt_id"),
  ("label_comp_id", "_atom_site.label_comp_id"),
  ("auth_asym_id", "_atom_site.auth_asym_id"),
  ("auth_seq_id", "_atom_site.auth_seq_id"),
  ("ins_code", "_atom_site.pdbx_PDB_ins_code"),
  ("cartn_x", "_atom_site.Cartn_x"),
  ("cartn_y", "_atom_site.Cartn_y"),
  ("cartn_z", "_atom_site.Cartn_z"),
  ("occupancy", "_atom_site.occupancy"),
  ("b_iso", "_atom_site.B_iso_or_equiv"),
  ("type_symbol", "_atom_site.type_symbol"),
  ("formal_charge", "_atom_site.pdbx_formal_charge"),
  ("resolution", "_atom_site.phenix_resolution"),
  ("fp", "_atom_site.phenix_scat_dispersion_real"),
  ("fdp", "_atom_site.phenix_scat_dispersion_imag"),
  ("label_asym_id", "_atom_site.label_asym_id"),
  ("label_entity_id", "_atom_site.label_entity_id"),
  ("label_seq_id", "_atom_site.label_seq_id"),
  ("auth_atom_id", "_atom_site.auth_atom_id"),
  ("model_num", "_atom_site.pdbx_PDB_model_num"),
  ("auth_segid", "_atom_site.auth_segid"),
  ("auth_break", "_atom_site.auth_break"),
)
aniso_tags = (
  ("aniso_id", "_atom_site_anisotrop.id"),
  ("aniso_auth_atom_id", "_atom_site_anisotrop.pdbx_auth_atom_id"),
  ("aniso_label_alt_id", "_atom_site_anisotrop.pdbx_label_alt_id"),
  ("aniso_auth_comp_id", "_atom_site_anisotrop.pdbx_auth_comp_id"),
  ("aniso_auth_asym_id", "_atom_site_anisotrop.pdbx_auth_asym_id"),
  ("aniso_auth_seq_id", "_atom_site_anisotrop.pdbx_auth_seq_id"),
  ("aniso_ins_code", "_atom_site_anisotrop.pdbx_PDB_ins_code"),
  ("aniso_u11", "_atom_site_anisotrop.U[1][1]"),
  ("aniso_u22", "_atom_site_anisotrop.U[2][2]"),
  ("aniso_u33", "_atom_site_anisotrop.U[3][3]"),
  ("aniso_u12", "_atom_site_anisotrop.U[1][2]"),
  ("aniso_u13", "_atom_site_anisotrop.U[1][3]"),
  ("aniso_u23", "_atom_site_anisotrop.U[2][3]"),
)
precision_names = ("coordinate_precision", "occupancy_precision",
  "b_iso_precision", "u_aniso_precision", "resolution_precision")
default_precisions = (5, 3, 5, 5, 2)


def cpp_labels(h, segid_as_auth_segid=False):
  """The arrays of h._atom_site_labels as an iotbx.pdb.atom_site_labels."""
  values = h._atom_site_labels(segid_as_auth_segid=segid_as_auth_segid)
  labels = iotbx.pdb.atom_site_labels()
  labels.chain_auth_asym_id = values.chain_auth_asym_id
  labels.rg_label_asym_id = values.rg_label_asym_id
  labels.rg_auth_seq_id = values.rg_auth_seq_id
  labels.ag_label_seq_id = values.ag_label_seq_id
  return labels


def cpp_columns(h, segid_as_auth_segid=False, output_break_records=False,
    **precisions):
  """atom_site_columns_from_hierarchy, with output_break_records resolved
  as as_cif_block resolves it."""
  args = [precisions.get(name, default)
    for name, default in zip(precision_names, default_precisions)]
  return iotbx.pdb.atom_site_columns_from_hierarchy(h,
    cpp_labels(h, segid_as_auth_segid=segid_as_auth_segid), *args,
    segid_as_auth_segid=segid_as_auth_segid,
    output_break_records=output_break_records and h.contains_break_records())


def check_block(h, cs=None, **kw):
  """Assert that as_cif_block writes the columns of
  atom_site_columns_from_hierarchy under their tags, in the order of
  atom_site_tags and aniso_tags: fp and fdp unless all ".", resolution
  unless all "0.00", auth_segid when asked for, auth_break when asked for
  and the hierarchy has break records, the anisotrop loop only when an atom
  has U, and the _chem_comp and _struct_asym ids. Returns the block."""
  out = cpp_columns(h, **kw)
  block = h.as_cif_block(cs, **kw)
  loop = block.get_loop("_atom_site")
  expected_tags = []
  for member, tag in atom_site_tags:
    column = list(getattr(out, member))
    if member in ("fp", "fdp") and column == ["."] * len(column):
      continue
    if member == "resolution" and column == ["0.00"] * len(column):
      continue
    if member == "auth_segid" and not kw.get("segid_as_auth_segid"):
      continue
    if member == "auth_break" and not (
        kw.get("output_break_records") and h.contains_break_records()):
      continue
    expected_tags.append(tag)
    assert list(loop[tag]) == column, tag
  assert list(loop.keys()) == expected_tags
  aniso = block.get_loop("_atom_site_anisotrop")
  if out.aniso_id.size() == 0:
    assert aniso is None
  else:
    assert list(aniso.keys()) == [tag for member, tag in aniso_tags]
    for member, tag in aniso_tags:
      assert list(aniso[tag]) == list(getattr(out, member)), tag
  for member, tag in (("chem_comp_ids", "_chem_comp.id"),
                      ("struct_asym_ids", "_struct_asym.id")):
    assert list(block.get(tag, [])) == list(getattr(out, member)), tag
  return block


def compare(name, h, cs, **kw):
  """check_block, and the same block as the Python reference."""
  block = check_block(h, cs, **kw)
  # PARITY-BEGIN
  from libtbx.test_utils import show_diff
  expected = str(h._as_cif_block_python(cs, **kw))
  assert not show_diff(str(block), expected), (name, kw)
  # PARITY-END


def compare_columns(name, h, cs=None, **kw):
  """The C++ columns of h. Each equals the column of the Python reference's
  block; a column the reference drops is all "." (fp, fdp), all "0.00"
  (resolution) or empty (auth_segid and auth_break when not written)."""
  out = cpp_columns(h, **kw)
  # PARITY-BEGIN
  block = h._as_cif_block_python(cs, **kw)
  loop = block.get_loop("_atom_site")
  n = loop.n_rows()
  for member, tag in atom_site_tags:
    column = list(getattr(out, member))
    if tag in loop:
      assert column == list(loop[tag]), (name, kw, tag, column)
    elif member in ("fp", "fdp"):
      assert column == ["."] * n, (name, kw, tag, column)
    elif member == "resolution":
      assert column == ["0.00"] * n, (name, kw, tag, column)
    else:
      assert column == [], (name, kw, tag, column)
  aniso = block.get_loop("_atom_site_anisotrop")
  for member, tag in aniso_tags:
    expected = [] if aniso is None else list(aniso[tag])
    assert list(getattr(out, member)) == expected, (name, kw, tag)
  for member, tag in (("chem_comp_ids", "_chem_comp.id"),
                      ("struct_asym_ids", "_struct_asym.id")):
    expected = list(block.get(tag, []))
    assert list(getattr(out, member)) == expected, (name, kw, tag)
  # PARITY-END
  return out


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
  out = cpp_columns(h)
  assert list(out.auth_asym_id)[-3:] == ["SEGA", "SEGA", "C"]
  assert list(out.label_asym_id)[-3:] == ["E", "E", "F"]


def exercise_cpp_columns(fixture_list):
  # every precision, each with its own value
  h, cs = fixture(fixture_list, "multi_model_pdb")
  out = cpp_columns(h, coordinate_precision=3, occupancy_precision=2,
    b_iso_precision=1, u_aniso_precision=4, resolution_precision=6)
  assert list(out.cartn_x[:2]) == ["1.000", "2.000"]
  assert list(out.occupancy[:1]) == ["1.00"]
  assert list(out.b_iso[:1]) == ["10.0"]
  assert list(out.aniso_u11) == ["0.0100", "0.0110"]
  assert list(out.resolution[:1]) == ["0.000000"]
  # _struct_asym in order of first appearance, not sorted
  h, cs = fixture(fixture_list, "many_ligands_pdb")
  out = cpp_columns(h)
  first_appearance = []
  for label_asym_id in out.label_asym_id:
    if label_asym_id not in first_appearance:
      first_appearance.append(label_asym_id)
  assert len(first_appearance) == 27
  assert list(out.struct_asym_ids) == first_appearance
  assert first_appearance != sorted(first_appearance)


def exercise_programmatic_hierarchy():
  # the default serial is "", which hy36decode rejects: one model raises
  root = build_root([""])
  try:
    cpp_columns(root)
  except ValueError:
    pass
  else:
    raise Exception_expected
  # PARITY-BEGIN
  try:
    root._as_cif_block_python()
  except ValueError:
    pass
  else:
    raise Exception_expected
  # PARITY-END
  # five blanks decode as 0; an empty model id is written as 1, a blank one
  # as the empty string; an empty insertion code as ?
  out = compare_columns("blank serial", build_root([""], serial="     "))
  assert list(out.id) == ["0", "0"]
  assert list(out.model_num) == ["1", "1"]
  assert list(out.type_symbol) == ["", ""]
  assert list(out.ins_code) == ["?", "?"]
  out = compare_columns("blank model id", build_root([" "], serial="     "))
  assert list(out.model_num) == ["", ""]
  # a hybrid-36 serial is decoded
  out = compare_columns("hybrid-36 serial", build_root([""], serial="A0000"))
  assert list(out.id) == ["100000", "100000"]
  # with two models the id is the row number; the serial is not read
  out = compare_columns("two models", build_root(["", "2"]))
  assert list(out.id) == ["1", "2", "3", "4"]
  assert list(out.model_num) == ["1", "1", "2", "2"]


def exercise_multi_model(fixture_list):
  h, cs = fixture(fixture_list, "multi_model_pdb")
  out = cpp_columns(h, segid_as_auth_segid=True, output_break_records=True)
  # serials restart in model 2: the id is the row, the anisotrop id too
  assert list(out.id) == [str(i) for i in range(1, 11)]
  assert list(out.aniso_id) == ["1", "9"]
  assert list(out.model_num) == ["1"] * 8 + ["2"] * 2
  assert list(out.auth_segid) == ["SEGA"] * 7 + ["    "] + ["SEGA"] * 2
  assert list(out.label_alt_id) == [".", ".", ".", "A", "B"] + ["."] * 5
  assert list(out.ins_code) == ["?", "?", "?", "A", "A"] + ["?"] * 5
  assert list(out.group_pdb) == ["ATOM"] * 6 + ["HETATM"] * 2 + ["ATOM"] * 2
  assert list(out.formal_charge) == ["?"] * 5 + ["1", "2", "?", "?", "?"]
  assert list(out.chem_comp_ids) == ["ALA", "GLY", "HOH", "SER", "ZN"]
  assert list(out.struct_asym_ids) == ["A", "B", "C", "E"]
  # the first atom after the BREAK and the first of chain B get 1; GLY 2A,
  # linked to ALA 1, and model 2, which starts the flag again, get 0
  assert list(out.auth_break) == ["0"] * 5 + ["1", "0", "1", "0", "0"]
  assert list(out.auth_atom_id) == [
    "N", "CA", "C", "N", "N", "N", "ZN", "O", "N", "CA"]
  assert list(out.label_entity_id) == ["?"] * 10
  assert list(out.label_asym_id) == ["A"] * 6 + ["B", "C", "E", "E"]
  assert list(out.auth_asym_id) == ["A"] * 7 + ["B", "A", "A"]
  assert list(out.auth_seq_id) == [
    "1", "1", "1", "2", "2", "5", "101", "1", "1", "1"]
  assert list(out.label_seq_id) == [
    "1", "1", "1", "2", "2", "3", ".", ".", "1", "1"]
  # an anisotrop row repeats the labels of its atom: U on the altloc A atom
  # of GLY 2A, which has an insertion code, and on the water of chain B
  h = h.deep_copy()
  h.atoms()[3].set_uij((0.1, 0.2, 0.3, 0.01, 0.02, 0.03))
  h.atoms()[7].set_uij((0.4, 0.5, 0.6, 0.04, 0.05, 0.06))
  out = compare_columns("multi_model_pdb, more U", h)
  assert list(out.aniso_id) == ["1", "4", "8", "9"]
  assert list(out.aniso_auth_atom_id) == ["N", "N", "O", "N"]
  assert list(out.aniso_label_alt_id) == [".", "A", ".", "."]
  assert list(out.aniso_auth_comp_id) == ["ALA", "GLY", "HOH", "ALA"]
  assert list(out.aniso_auth_asym_id) == ["SEGA", "SEGA", "B", "SEGA"]
  assert list(out.aniso_auth_seq_id) == ["1", "2", "1", "1"]
  assert list(out.aniso_ins_code) == ["?", "A", "?", "?"]
  assert list(out.aniso_u11) == ["0.01000", "0.10000", "0.40000", "0.01100"]
  assert list(out.aniso_u22) == ["0.02000", "0.20000", "0.50000", "0.02100"]
  assert list(out.aniso_u33) == ["0.03000", "0.30000", "0.60000", "0.03100"]
  assert list(out.aniso_u12) == ["0.00100", "0.01000", "0.04000", "0.00110"]
  assert list(out.aniso_u13) == ["0.00200", "0.02000", "0.05000", "0.00210"]
  assert list(out.aniso_u23) == ["0.00300", "0.03000", "0.06000", "0.00310"]


chain_id_fallbacks_pdb = """\
ATOM      1  N   ALA     1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  N   GLY B   1       4.000   2.000   3.000  1.00 10.00      LONG N
ATOM      3  N   SER C   1       7.000   2.000   3.000  1.00 10.00      LO   N
END
"""


def exercise_chain_id_fallbacks():
  inp = iotbx.pdb.input(source_info=None,
    lines=flex.split_lines(chain_id_fallbacks_pdb))
  h = inp.construct_hierarchy()
  # a blank chain id with a blank segid falls back to the label_asym_id; a
  # segid longer than the chain id is the auth_asym_id unless it is written
  # as auth_segid
  out = compare_columns("chain id fallbacks", h)
  assert list(out.auth_asym_id) == ["A", "LONG", "LO"]
  out = compare_columns("chain id fallbacks", h, segid_as_auth_segid=True)
  assert list(out.auth_asym_id) == ["A", "B", "C"]
  assert list(out.auth_segid) == ["    ", "LONG", "LO  "]


def charge_root():
  """Seven atoms of one residue, one per formal charge spelling."""
  root = build_root([""], serial="    1")
  ag = root.only_atom_group()
  for atom in ag.atoms():
    ag.remove_atom(atom)
  for i, charge in enumerate(["2+", "+2", "1-", "-1", "0", "  ", "?"]):
    atom = iotbx.pdb.hierarchy.atom()
    atom.set_name(" C%d " % i)
    atom.serial = "%5d" % (i + 1)
    atom.set_charge(charge)
    ag.append_atom(atom)
  return root


def exercise_charge_fp_resolution():
  root = charge_root()
  # charge_tidy, then the sign in front; fp and fdp on one atom keep both
  # columns; resolution all zero drops its column
  atom = root.atoms()[2]
  atom.set_fp(-0.5)
  atom.set_fdp(1.25)
  out = compare_columns("charges", root)
  assert list(out.formal_charge) == ["2", "2", "-1", "-1", "?", "?", "?"]
  assert list(out.fp) == [".", ".", "-0.5000", ".", ".", ".", "."]
  assert list(out.fdp) == [".", ".", "1.2500", ".", ".", ".", "."]
  loop = root.as_cif_block().get_loop("_atom_site")
  assert "_atom_site.phenix_scat_dispersion_real" in loop
  assert "_atom_site.phenix_resolution" not in loop
  # fdp alone; a resolution on one atom keeps the column
  root.atoms()[4].set_fdp(0.3)
  root.atoms()[5].set_resolution(2.5)
  out = compare_columns("charges", root)
  assert list(out.fp)[4] == "0.0000"
  assert list(out.resolution)[4:6] == ["0.00", "2.50"]
  # with atoms, a resolution precision other than 2 never matches the
  # "0.00" test
  out = compare_columns("charges", charge_root(), resolution_precision=3)
  assert list(out.resolution) == ["0.000"] * 7


break_records_pdb = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00 10.00           C
BREAK
ATOM      3  N   GLY A   2       5.000   2.000   3.000  1.00 10.00           N
ATOM      4  CA  GLY A   2       6.000   2.000   3.000  1.00 10.00           C
TER
ATOM      5  N   SER B   1       9.000   2.000   3.000  1.00 10.00           N
END
"""


def exercise_break_records():
  inp = iotbx.pdb.input(source_info=None,
    lines=flex.split_lines(break_records_pdb))
  h = inp.construct_hierarchy()
  assert [rg.link_to_previous for rg in h.residue_groups()] == [
    False, False, False]
  # the first atom of a residue group after a break gets 1; the first
  # residue group of chain B does too, because the first-in-chain flag of
  # the Python reference is reset per model, not per chain
  out = compare_columns("break records", h, output_break_records=True)
  assert list(out.auth_break) == ["0", "0", "1", "0", "1"]
  out = compare_columns("break records", h)
  assert list(out.auth_break) == []


def exercise_python_whitespace():
  # str.strip() removes the Unicode blanks too; the hierarchy holds them as
  # UTF-8
  root = build_root(["\u2003 2"], serial="    1")
  atom = root.atoms()[0]
  atom.set_name("N\xa0")
  atom.set_element("\x1fN")
  root.only_atom_group().resname = "\u3000ALA\x1c"
  out = compare_columns("python whitespace", root)
  assert list(out.label_atom_id) == ["N", "CA"]
  assert list(out.type_symbol) == ["N", ""]
  assert list(out.label_comp_id) == ["ALA", "ALA"]
  assert list(out.model_num) == ["2", "2"]
  # every str.isspace() character, at both ends of the model id, the text
  # field without a length limit
  blanks = ""
  for i in range(0x3001):
    if chr(i).isspace():
      blanks += chr(i)
  assert len(blanks) == 29
  out = compare_columns("all blanks",
    build_root([blanks + "7" + blanks], serial="    1"))
  assert list(out.model_num) == ["7", "7"]


def exercise_non_finite_and_decimal_point():
  root = build_root([""], serial="    1")
  atoms = root.atoms()
  atoms[0].set_xyz((float("nan"), -float("nan"), float("inf")))
  atoms[1].set_xyz((-float("inf"), -0.0, 0.125))
  out = compare_columns("non-finite", root, coordinate_precision=2)
  assert list(out.cartn_x) == ["nan", "-inf"]
  assert list(out.cartn_y) == ["nan", "-0.00"]
  assert list(out.cartn_z) == ["inf", "0.12"]
  # the C library writes the decimal point of LC_NUMERIC; Python's % and
  # the writer always write '.', and the writer leaves the locale alone
  previous = locale.setlocale(locale.LC_NUMERIC)
  comma_locale = None
  for name in ("de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "fr_FR.utf8",
               "en_DK.UTF-8", "en_DK.utf8"):
    try:
      locale.setlocale(locale.LC_NUMERIC, name)
    except locale.Error:
      continue
    comma_locale = name
    break
  if comma_locale is None:
    print("exercise_non_finite_and_decimal_point: no decimal comma locale")
    return
  try:
    assert locale.localeconv()["decimal_point"] == ","
    out = compare_columns("decimal comma", root, coordinate_precision=2)
    assert list(out.cartn_z) == ["inf", "0.12"]
    root.as_cif_block()
    assert locale.localeconv()["decimal_point"] == ","
  finally:
    locale.setlocale(locale.LC_NUMERIC, previous)


def exercise_rare_branches():
  # an atom group without atoms beside one with atoms: no _chem_comp id
  # for it; a 70-digit coordinate, longer than the formatting buffer; a
  # two-byte blank leading a name and a three-byte blank ending a residue
  # name
  root = build_root([""], serial="    1")
  root.only_residue_group().append_atom_group(
    iotbx.pdb.hierarchy.atom_group(altloc="B", resname="GLY"))
  atom = root.atoms()[0]
  atom.set_xyz((1e70, -1e70, 0.5))
  atom.set_name("\xa0N")
  list(root.atom_groups())[0].resname = "ALA\u3000"
  compare("rare branches", root, None)
  out = cpp_columns(root)
  assert list(out.chem_comp_ids) == ["ALA"]
  assert len(out.cartn_x[0]) == 77 and out.cartn_x[0].endswith(".00000")
  assert list(out.label_atom_id) == ["N", "CA"]
  assert list(out.label_comp_id) == ["ALA", "ALA"]
  # a residue group whose only atom group has no atoms: its label_asym_id
  # is not written to _struct_asym
  root = build_root([""], serial="    1")
  rg = iotbx.pdb.hierarchy.residue_group(resseq="   2")
  rg.append_atom_group(iotbx.pdb.hierarchy.atom_group(resname="HOH"))
  root.only_chain().append_residue_group(rg)
  out = compare_columns("empty residue group", root)
  assert list(cpp_labels(root).rg_label_asym_id) == ["A", "B"]
  assert list(out.struct_asym_ids) == ["A"]


def exercise_errors(fixture_list):
  h, cs = fixture(fixture_list, "multi_model_pdb")
  args = default_precisions + (False, False)
  # a labels array of the wrong length is named
  for member in ("chain_auth_asym_id", "rg_label_asym_id", "rg_auth_seq_id",
                 "ag_label_seq_id"):
    for change in ("short", "long"):
      labels = cpp_labels(h)
      values = getattr(labels, member)
      if change == "short":
        values = values[:-1]
      else:
        values.append("X")
      setattr(labels, member, values)
      try:
        iotbx.pdb.atom_site_columns_from_hierarchy(h, labels, *args)
      except ValueError as e:
        assert member in str(e), str(e)
      else:
        raise Exception_expected
  # the documented ValueError: a negative precision, as Python's % raises,
  # once a value is formatted with it; a serial or a four-character residue
  # number that is not hybrid-36
  out = cpp_columns(iotbx.pdb.hierarchy.root(), coordinate_precision=-1)
  assert out.cartn_x.size() == 0
  bad_resseq = build_root([""], serial="    1")
  bad_resseq.only_residue_group().resseq = "1A2B"
  for writer in (lambda: h.as_cif_block(coordinate_precision=-1),
                 lambda: build_root([""], serial="  x  ").as_cif_block(),
                 bad_resseq.as_cif_block):
    try:
      writer()
    except ValueError:
      pass
    else:
      raise Exception_expected
  # PARITY-BEGIN
  try:
    h._as_cif_block_python(coordinate_precision=-1)
  except ValueError:
    pass
  else:
    raise Exception_expected
  # PARITY-END
  # the documented IndexError, and the difference in order: the labels of
  # all chains are resolved before any atom is written, so the empty chain
  # B raises before the bad serial of chain A, which the Python body hit
  # first
  root = build_root([""], serial="  x  ")
  root.only_model().append_chain(iotbx.pdb.hierarchy.chain(id="B"))
  try:
    root.as_cif_block()
  except IndexError:
    pass
  else:
    raise Exception_expected
  # PARITY-BEGIN
  try:
    root._as_cif_block_python()
  except ValueError:
    pass
  else:
    raise Exception_expected
  # PARITY-END


def exercise_input_unchanged(fixture_list):
  # the writer reads its input only: i_seq, tmp, serials, names and the
  # atom order stay as they were (the copy's i_seq are all 0)
  h, cs = fixture(fixture_list, "multi_model_pdb")
  h = h.select(flex.bool(h.atoms_size(), True), copy_atoms=True)
  atoms = h.atoms()
  for i, atom in enumerate(atoms):
    atom.tmp = 100 + i
  def state():
    atoms = h.atoms()
    return (list(atoms.extract_i_seq()), list(atoms.extract_tmp_as_size_t()),
      list(atoms.extract_serial()), list(atoms.extract_name()),
      list(atoms.extract_xyz()))
  before = state()
  assert before[0] == [0] * 10
  h.as_cif_block(cs, segid_as_auth_segid=True, output_break_records=True)
  assert state() == before


def exercise_atoms_with_other_parent():
  # The C++ writer takes label_alt_id from the atom group that holds the
  # atom; the Python reference took it from atom.parent(). An atom appended
  # with append_atom_with_other_parent keeps its first group as parent
  h = iotbx.pdb.input(source_info=None,
    lines=flex.split_lines(multi_model_pdb)).construct_hierarchy()
  gly_a = list(h.atom_groups())[1]
  z = iotbx.pdb.hierarchy.atom_group(altloc="Z", resname="GLY")
  z.append_atom_with_other_parent(gly_a.atoms()[0])
  rg = iotbx.pdb.hierarchy.residue_group(resseq="   9")
  rg.append_atom_group(z)
  h.models()[1].chains()[0].append_residue_group(rg)
  assert cpp_columns(h).label_alt_id[-1] == "Z"
  # PARITY-BEGIN
  reference = h._as_cif_block_python().get_loop("_atom_site")
  assert reference["_atom_site.label_alt_id"][-1] == "A"
  # PARITY-END
  # select() without copy_atoms keeps each atom's parent, the atom group of
  # the hierarchy selected from: the same value while the original is
  # unchanged
  h = iotbx.pdb.input(source_info=None,
    lines=flex.split_lines(multi_model_pdb)).construct_hierarchy()
  selected = h.select(flex.bool(h.atoms_size(), True))
  compare("selection", selected, None)
  # an altloc changed after the select is written from the selection (the
  # reference wrote the original's)
  ag = list(selected.atom_groups())[1]
  assert ag.altloc == "A"
  ag.altloc = "C"
  assert cpp_columns(selected).label_alt_id[3] == "C"
  # PARITY-BEGIN
  reference = selected._as_cif_block_python().get_loop("_atom_site")
  assert reference["_atom_site.label_alt_id"][3] == "A"
  # PARITY-END
  ag.altloc = "A"
  # once the original is gone the parent is None: the block is written
  # (the reference raised AttributeError)
  del h
  assert selected.atoms()[0].parent() is None
  # PARITY-BEGIN
  try:
    selected._as_cif_block_python()
  except AttributeError:
    pass
  else:
    raise Exception_expected
  # PARITY-END
  out = cpp_columns(selected)
  assert list(out.label_alt_id) == [".", ".", ".", "A", "B"] + ["."] * 5
  check_block(selected)


def exercise_switched(fixture_list):
  # PARITY-BEGIN
  # as_cif_block does not go through the Python reference
  reference = iotbx.pdb.hierarchy.root._as_cif_block_python
  def refuse(self, *args, **kw):
    raise AssertionError("as_cif_block called _as_cif_block_python")
  h, cs = fixture(fixture_list, "multi_model_pdb")
  iotbx.pdb.hierarchy.root._as_cif_block_python = refuse
  try:
    h.as_cif_block()
  finally:
    iotbx.pdb.hierarchy.root._as_cif_block_python = reference
  # PARITY-END
  # every fixture with the defaults; multi_model_pdb, which has segids and
  # break records, with every other combination of the flags and of a
  # coordinate precision other than the default
  for name, h, cs in fixture_list:
    compare(name, h, cs)
  h, cs = fixture(fixture_list, "multi_model_pdb")
  for segid in (False, True):
    for breaks in (False, True):
      for precision in (5, 3):
        if (segid, breaks, precision) != (False, False, 5):
          compare("multi_model_pdb", h, cs, segid_as_auth_segid=segid,
            output_break_records=breaks, coordinate_precision=precision)
  # an auth_segid column of blanks, fp and fdp with a coordinate precision
  # other than the default; no break records: no auth_break column even
  # when asked for
  h, cs = fixture(fixture_list, "input_anom")
  compare("input_anom", h, cs, segid_as_auth_segid=True,
    coordinate_precision=3)
  h, cs = fixture(fixture_list, "input_1ab1")
  assert not h.contains_break_records()
  compare("input_1ab1", h, cs, output_break_records=True)
  # every precision through as_cif_block, each with its own value, and a
  # kept resolution column
  h, cs = fixture(fixture_list, "multi_model_pdb")
  compare("multi_model_pdb", h, None, coordinate_precision=3,
    occupancy_precision=2, b_iso_precision=1, u_aniso_precision=4,
    resolution_precision=6)
  root = charge_root()
  root.atoms()[5].set_resolution(2.5)
  compare("resolution", root, None)


def exercise_block_layout(fixture_list):
  # the loops in order, the symmetry items, the _chem_comp and _struct_asym
  # ids; no symmetry items without a crystal_symmetry
  h, cs = fixture(fixture_list, "multi_model_pdb")
  block = h.as_cif_block(cs)
  assert list(block.loops.keys()) == ["_atom_site", "_atom_site_anisotrop",
    "_space_group_symop", "_chem_comp.id", "_struct_asym.id"]
  assert block["_cell.length_a"] == "20.000"
  assert block["_space_group.name_H-M_alt"] == "P 1"
  assert list(block["_chem_comp.id"]) == ["ALA", "GLY", "HOH", "SER", "ZN"]
  assert list(block["_struct_asym.id"]) == ["A", "B", "C", "E"]
  block = h.as_cif_block()
  assert "_cell.length_a" not in block
  assert "_space_group_symop" not in block.loops


def run():
  # all_label_asym_ids() builds its list of 475254 ids, A to ZZZZ, on every
  # call, 0.2 s, and every as_cif_block call makes one call; the checks here
  # make hundreds, so they share one list, up to ZZZ, more than any of them
  # needs.
  all_label_asym_ids = iotbx.pdb.hierarchy.all_label_asym_ids
  label_asym_ids = all_label_asym_ids(maximum_length=3)
  iotbx.pdb.hierarchy.all_label_asym_ids = lambda: label_asym_ids
  try:
    fixture_list = fixtures()
    exercise_atom_site_labels(fixture_list)
    exercise_cpp_columns(fixture_list)
    exercise_programmatic_hierarchy()
    exercise_multi_model(fixture_list)
    exercise_chain_id_fallbacks()
    exercise_charge_fp_resolution()
    exercise_break_records()
    exercise_python_whitespace()
    exercise_non_finite_and_decimal_point()
    exercise_rare_branches()
    exercise_errors(fixture_list)
    exercise_input_unchanged(fixture_list)
    exercise_atoms_with_other_parent()
    exercise_switched(fixture_list)
    exercise_block_layout(fixture_list)
  finally:
    iotbx.pdb.hierarchy.all_label_asym_ids = all_label_asym_ids


if __name__ == "__main__":
  run()
  print("OK")
