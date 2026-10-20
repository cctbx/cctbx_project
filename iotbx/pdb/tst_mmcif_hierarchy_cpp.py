from __future__ import absolute_import, division, print_function
"""
The C++ mmCIF hierarchy builder (iotbx.pdb.build_mmcif_hierarchy on the columns
of iotbx.pdb.mmcif.atom_site_columns_from_block) against the Python builder
(iotbx.pdb.mmcif.pdb_hierarchy_builder_python): the column binding, cell
parsing, the boundary, naming and field rules, and each difference F1-F11
listed in iotbx/pdb/mmcif_hierarchy.h, asserted on both builders.
"""
import sys
import math
from six.moves import cStringIO as StringIO
from libtbx.test_utils import show_diff
import iotbx.cif
from iotbx.pdb import mmcif
from iotbx.pdb import atom_site_columns, build_mmcif_hierarchy
from iotbx.pdb import mmcif_parse_cell_double, mmcif_parse_cell_int
from iotbx.pdb.tst_mmcif_hierarchy import input_1ab1, input_4edr, input_null_occ_b
from iotbx.pdb.tst_mmcif_hierarchy_parity import compare_hierarchies
from cctbx.array_family import flex
from cctbx import adptbx


def exercise_cpp_columns_skeleton():
  cols = atom_site_columns()
  assert cols.type_symbol.size() == 0
  assert cols.aniso_is_b is False
  cols.type_symbol = flex.std_string(["N"])
  assert list(cols.type_symbol) == ["N"]
  # zero rows: an empty hierarchy, no exception
  h = build_mmcif_hierarchy(atom_site_columns())
  assert h.models_size() == 0
  # a mandatory column of the wrong length is reported by name
  cols = atom_site_columns()
  cols.atom_name = flex.std_string(["N", "CA"])
  cols.type_symbol = flex.std_string(["N"])
  try:
    build_mmcif_hierarchy(cols)
  except ValueError as e:
    assert "type_symbol" in str(e) and "2" in str(e), str(e)
  else:
    raise AssertionError("length mismatch not detected")
  # a present optional column of the wrong length is reported too
  cols = atom_site_columns()
  for name in ("type_symbol", "atom_name", "alt_id", "label_asym_id",
               "auth_asym_id", "comp_id", "seq_id", "cartn_x", "cartn_y",
               "cartn_z", "b_iso", "occupancy"):
    setattr(cols, name, flex.std_string(["x"]))
  cols.segid = flex.std_string(["A", "B"])
  try:
    build_mmcif_hierarchy(cols)
  except ValueError as e:
    assert "segid" in str(e), str(e)
  else:
    raise AssertionError("optional length mismatch not detected")

def exercise_cpp_cell_parsing():
  # everything flex.double accepts, with the same value bit for bit
  # (no subnormal literal here: libstdc++'s num_get accepts 4.9e-324 while
  # libc++'s rejects it, so underflow is unspecified and not tested; overflow
  # is rejected everywhere and is in the reject list below)
  for cell in ("1.5", "-2", "+3.25", ".5", "5.", "1e3", "1E-3", "12.345(6)",
               "1.5+03", "1.5-03", "0", "-0.0", "inf", "nan", "0x1p-2",
               "1.7976931348623157e308"):
    try:
      expected = flex.double(flex.std_string([cell]))[0]
    except Exception:
      continue   # a cell flex.double rejects is skipped here
    got = mmcif_parse_cell_double(cell)
    assert repr(got) == repr(expected), (cell, got, expected)
  # everything flex.double rejects
  for cell in ("", " 1.5", "1.5 ", "1.5x", "?", ".", "1e999", "1.5(", "(5)",
               "1.5()", "1.5(a)", "abc", "1,5", "1.5d0", "1.5(99999999999)",
               "1.5(2147483647)", "1.5(2147483648)", "1.5(+5)", "1.5(-5)",
               "1.5\x00abc"):
    try:
      flex.double(flex.std_string([cell]))
      flex_ok = True
    except Exception:
      flex_ok = False
    try:
      mmcif_parse_cell_double(cell)
      ours_ok = True
    except ValueError:
      ours_ok = False
    assert flex_ok == ours_ok, (cell, flex_ok, ours_ok)
  assert mmcif_parse_cell_int("12") == 12
  assert mmcif_parse_cell_int(" -7 ") == -7
  assert mmcif_parse_cell_int("+5") == 5
  for cell in ("", "1.5", "1a", "a1", "?", ".", "1 2", "--1", "99999999999999999999"):
    try:
      mmcif_parse_cell_int(cell)
    except ValueError:
      pass
    else:
      raise AssertionError("accepted %r" % cell)

atom_site_header = """\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""
# column order of the rows below:
# group id elem name alt comp lasym ent lseq icode x y z occ b aseq acomp aasym aname model


def _block(text):
  return list(iotbx.cif.reader(input_string=text).model().values())[0]


def _drop_column(text, tag):
  """Remove tag and its field from every ATOM/HETATM row of a fixture."""
  lines = text.splitlines()
  tags = [l for l in lines if l.startswith("_atom_site.")]
  index = tags.index(tag)
  out = []
  for l in lines:
    if l == tag:
      continue
    if l.startswith("ATOM") or l.startswith("HETATM"):
      f = l.split()
      del f[index]
      l = " ".join(f)
    out.append(l)
  return "\n".join(out) + "\n"


def _python(block):
  """The Python builder, the reference."""
  return mmcif.pdb_hierarchy_builder_python(block).hierarchy


def _cpp_block(block):
  """The C++ builder called directly on the columns of
  atom_site_columns_from_block."""
  return build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(block))


def _cpp(text):
  return _cpp_block(_block(text))


def _both(text):
  """(Python build, C++ build), both unsorted."""
  return (_python(_block(text)), _cpp(text))


def exercise_cpp_core_parity():
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 1
ATOM 2 C CA A ALA A 1 5 ? 1.5 2.0 3.0 0.50 11.00 5 ALA A CA 1
ATOM 3 C CA B ALA A 1 5 ? 1.6 2.0 3.0 0.50 12.00 5 ALA A CA 1
ATOM 4 N N . GLY A 1 6 ? 2.0 2.0 3.0 1.00 13.00 6 GLY A N 1
ATOM 5 N N . GLY A 1 6 A 2.5 2.0 3.0 1.00 13.00 6 GLY A N 1
HETATM 6 O O . HOH B 2 . ? 9.0 9.0 9.0 1.00 30.00 101 HOH A O 1
HETATM 7 FE FE . HEM C 3 . ? 8.0 9.0 9.0 1.00 30.00 201 HEM A FE 1
ATOM 8 N N . SER D 4 1 ? 0.0 0.0 0.0 1.00 5.00 1 SER B N 1
ATOM 9 N N . SER D 4 1 ? 0.0 0.0 0.0 1.00 5.00 1 SER B N 2
"""
  block = _block(text)
  ref = _python(block)
  new = build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(block))
  compare_hierarchies(ref, new)
  assert new.models_size() == 2
  assert [c.id for c in new.models()[0].chains()] == ["A", "A", "B"]
  rgs = new.models()[0].chains()[0].residue_groups()
  assert [rg.resid() for rg in rgs] == ["   5 ", "   6 ", "   6A"]
  assert [ag.altloc for ag in rgs[0].atom_groups()] == ["", "A", "B"]
  assert new.models()[0].chains()[0].atoms()[0].name == " N  "
  assert new.models()[0].chains()[1].atoms()[1].name == "FE  "
  assert new.models()[0].chains()[1].atoms()[0].hetero
  # two blank-altloc residue names in one residue group: each new blank group
  # is inserted at the front, so the order is reversed (as in the Python builder)
  text = atom_site_header + """\
HETATM 1 O O . HOH A 1 . ? 1.0 2.0 3.0 1.00 10.00 7 HOH A O 1
HETATM 2 NA NA . NA A 1 . ? 1.0 2.0 3.0 1.00 10.00 7 NA A NA 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  rg0 = list(new.residue_groups())[0]   # root.residue_groups() is a generator
  assert [ag.resname for ag in rg0.atom_groups()] == [" NA", "HOH"]


def exercise_zero_rows():
  # _atom_site present with zero rows -> empty hierarchy, also when
  # _atom_site_anisotrop has rows (there is no atom to join them to)
  for text in (atom_site_header,
               atom_site_header + aniso_u_loop + "1 0.1 0.2 0.3 0.01 0.02 0.03\n"):
    ref, new = _both(text)
    assert ref.models_size() == 0 and new.models_size() == 0
    compare_hierarchies(ref, new)


def exercise_model_boundaries():
  # a new model starts wherever the model number changes, so non-adjacent
  # repeats of a number are separate models
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 1
ATOM 2 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 2
ATOM 3 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [m.id for m in new.models()] == ["1", "2", "1"]
  # single model: id becomes ""
  ref, new = _both(atom_site_header + "ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 7\n")
  compare_hierarchies(ref, new)
  assert new.models()[0].id == ""


def exercise_chain_break_rule():
  """
  The chain rules: a new chain where the model or auth_asym_id changes, or
  where label_asym_id changes after a polymer residue (rule 3).
  """
  # polymer, then ligand and water in the same auth chain: split after the polymer only
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
ATOM 2 N N . GLY A 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1
HETATM 3 C C1 . LIG B 2 . ? 1.0 2.0 3.0 1.00 10.00 101 LIG A C1 1
HETATM 4 C C1 . LIG C 2 . ? 1.0 2.0 3.0 1.00 10.00 102 LIG A C1 1
HETATM 5 O O . HOH D 3 . ? 1.0 2.0 3.0 1.00 10.00 201 HOH A O 1
HETATM 6 O O . HOH D 3 . ? 1.0 2.0 3.0 1.00 10.00 202 HOH A O 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [len(c.residue_groups()) for c in new.models()[0].chains()] == [2, 4]
  # two same-id protein chains as the Phenix writer emits them: split
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
ATOM 2 N N . GLY A 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1
ATOM 3 N N . ALA B 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
ATOM 4 N N . GLY B 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [c.id for c in new.models()[0].chains()] == ["A", "A"]
  # a D-amino acid before the label change: no split today, reproduced
  text = atom_site_header + """\
ATOM 1 N N . DAL A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 DAL A N 1
ATOM 2 N N . ALA B 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 ALA A N 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert len(new.models()[0].chains()) == 1
  # auth_break clears link_to_previous on a later residue group only
  text = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.auth_break\n") + """\
ATOM 1 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1 1
ATOM 2 N N . GLY A 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1 0
ATOM 3 N N . SER A 1 3 ? 1.0 2.0 3.0 1.00 10.00 3 SER A N 1 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [rg.link_to_previous for rg in new.residue_groups()] == [True, True, False]


def exercise_atom_name_and_element():
  # a name shorter than four characters gets a leading space when the element
  # has one letter and the name starts with it (CA of carbon, not of calcium)
  text = atom_site_header + """\
ATOM 1 C CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1
ATOM 2 FE FE1 . HEM A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A FE1 1
ATOM 3 H HD21 . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A HD21 1
ATOM 4 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1
ATOM 5 CA CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [a.name for a in new.atoms()] == [" CA ", "FE1 ", "HD21", " C  ", "CA  "]
  # over-long fields: ValueError naming the tag and the row
  for bad_row, what, tag in (
      ("ATOM 1 C HD211 . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A HD211 1",
       "name", "label_atom_id"),
      ("ATOM 1 FEX FE . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A FE 1",
       "element", "_atom_site.type_symbol"),
      ("ATOM 1 C C . ALA A 1 1 AB 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1",
       "icode", "pdbx_PDB_ins_code")):
    try:
      _cpp(atom_site_header + bad_row + "\n")
    except ValueError as e:
      assert tag in str(e) and "row 1" in str(e), str(e)
    else:
      raise AssertionError("over-long %s accepted" % what)
  # a multi-character altloc is accepted by both builders (atom_group altloc
  # is a std::string, not a fixed-width field)
  ref, new = _both(atom_site_header +
    "ATOM 1 C C AB ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n")
  compare_hierarchies(ref, new)
  assert list(new.residue_groups())[0].atom_groups()[0].altloc == "AB"
  # F8: empty name after stripping -> four spaces (the Python builder raises IndexError)
  text = atom_site_header + "ATOM 1 C ' ' . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A ' ' 1\n"
  try:
    _python(_block(text))
  except IndexError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; update divergence F8 in iotbx/pdb/mmcif_hierarchy.h")
  assert _cpp(text).atoms()[0].name == "    "


def exercise_charge():
  # formal charge: the sign comes from a leading or trailing '-', zero has no
  # sign, ? and . leave the charge blank
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.pdbx_formal_charge\n")
  rows = ""
  for i, charge in enumerate(("2+", "+2", "-1", "1-", "0", "?", ".")):
    rows += "ATOM %d C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 %s\n" % (i + 1, charge)
  ref, new = _both(header + rows)
  compare_hierarchies(ref, new)
  assert [a.charge for a in new.atoms()] == ["2+", "2+", "1-", "1-", "0", "", ""]
  try:
    _cpp(header + "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 x+\n")
  except ValueError as e:
    assert "pdbx_formal_charge" in str(e) and "row 1" in str(e), str(e)
  else:
    raise AssertionError("bad charge accepted")
  try:
    _cpp(header + "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 10+\n")
  except ValueError as e:
    assert "pdbx_formal_charge" in str(e) and "row 1" in str(e), str(e)
  else:
    raise AssertionError("over-long charge accepted")


def exercise_segid_and_serial():
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.auth_segid\n")
  text = header + """\
ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 AB
ATOM 100000 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 ABCDEF
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [a.segid for a in new.atoms()] == ["AB  ", "ABCD"]
  assert [a.serial for a in new.atoms()] == ["    1", "A0000"]
  # resseq hybrid-36 and the four-character fallback
  text = atom_site_header + """\
ATOM 1 C C . ALA A 1 10000 ? 1.0 2.0 3.0 1.00 10.00 10000 ALA A C 1
ATOM 2 C C . ALA A 1 -5 ? 1.0 2.0 3.0 1.00 10.00 -5 ALA A C 1
ATOM 3 C C . ALA A 1 XYZQ ? 1.0 2.0 3.0 1.00 10.00 XYZQ ALA A C 1
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert [rg.resseq for rg in new.residue_groups()] == ["A000", "  -5", "XYZQ"]


def exercise_null_rejection_and_nulls_in_b_occ():
  # occupancy and B accept . and ? as 0.0 (integrative models write them);
  # coordinates do not
  text = atom_site_header + "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 ? . 1 ALA A C 1\n"
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  assert new.atoms()[0].occ == 0.0 and new.atoms()[0].b == 0.0
  for cell_row in ("ATOM 1 C C . ALA A 1 1 ? ? 2.0 3.0 1.0 1.0 1 ALA A C 1",
                   "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.0 1.0.1 1 ALA A C 1"):
    try:
      _cpp(atom_site_header + cell_row + "\n")
    except ValueError as e:
      assert "row 1" in str(e), str(e)
    else:
      raise AssertionError("bad numeric cell accepted")


aniso_u_loop = """\
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
"""


def exercise_aniso_join():
  rows = """\
ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1
ATOM 2 C CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1
ATOM 3 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
"""
  # ids 1 and 2 have U; id 2 twice (the last row wins, as in a Python dict);
  # id 3 none
  text = atom_site_header + rows + aniso_u_loop + """\
1 0.1 0.2 0.3 0.01 0.02 0.03
2 0.4 0.5 0.6 0.04 0.05 0.06
2 0.7 0.8 0.9 0.07 0.08 0.09
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  u = [a.uij for a in new.atoms()]
  assert u[0] == (0.1, 0.2, 0.3, 0.01, 0.02, 0.03)
  assert u[1] == (0.7, 0.8, 0.9, 0.07, 0.08, 0.09)
  assert u[2] == (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0)
  # F1: id above 99999 is joined by the id string; the Python builder loses it
  big = atom_site_header + "ATOM 123456 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n" + \
    aniso_u_loop + "123456 0.1 0.2 0.3 0.01 0.02 0.03\n"
  block = _block(big)
  assert _python(block).atoms()[0].uij == (-1.0,) * 6
  assert build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(block)).atoms()[0].uij == \
    (0.1, 0.2, 0.3, 0.01, 0.02, 0.03)
  # F2: B tensor converted to U (the Python builder raises NameError)
  textb = atom_site_header + rows[:rows.index("ATOM 2")] + aniso_u_loop.replace("U[", "B[") + \
    "1 8.0 8.0 8.0 0.0 0.0 0.0\n"
  block = _block(textb)
  try:
    _python(block)
  except NameError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; update divergence F2 in iotbx/pdb/mmcif_hierarchy.h")
  h = build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(block))
  expected = adptbx.b_as_u((8.0, 8.0, 8.0, 0.0, 0.0, 0.0))
  assert all(abs(x - y) < 1e-12 for x, y in zip(h.atoms()[0].uij, expected))
  # F4 with aniso: no _atom_site.id is an error
  text4 = _drop_column(atom_site_header + rows[:rows.index("ATOM 2")], "_atom_site.id") + \
    aniso_u_loop + "1 0.1 0.2 0.3 0.01 0.02 0.03\n"
  try:
    build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(_block(text4)))
  except ValueError as e:
    assert "_atom_site.id" in str(e), str(e)
  else:
    raise AssertionError("aniso without atom ids accepted")
  # a bad aniso cell names the anisotrop row
  textbad = atom_site_header + rows + aniso_u_loop + "1 0.1 0.2 x 0.01 0.02 0.03\n"
  try:
    build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(_block(textbad)))
  except ValueError as e:
    assert "_atom_site_anisotrop" in str(e) and "row 1" in str(e), str(e)
  else:
    raise AssertionError("bad aniso cell accepted")


def exercise_fp_fdp_resolution_cpp():
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.phenix_scat_dispersion_real\n"
    "_atom_site.phenix_scat_dispersion_imag\n_atom_site.phenix_resolution\n")
  text = header + """\
ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 0.5 1.5 2.0
ATOM 2 C CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1 ? . ?
"""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  a = new.atoms()
  assert (a[0].fp, a[0].fdp, a[0].resolution) == (0.5, 1.5, 2.0)
  assert (a[1].fp, a[1].fdp, a[1].resolution) == (0.0, 0.0, 0.0)


def exercise_divergences_f3_f5_f6_f7_f9():
  one = "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n"
  # F3: no model number column -> one model with id ""; Python: TypeError
  text = _drop_column(atom_site_header + one, "_atom_site.pdbx_PDB_model_num")
  try:
    _python(_block(text))
  except TypeError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; update divergence F3 in iotbx/pdb/mmcif_hierarchy.h")
  h = _cpp_block(_block(text))
  assert h.models_size() == 1 and h.models()[0].id == ""
  # F5: missing type_symbol or numeric column -> ValueError naming the tag
  for tag in ("_atom_site.type_symbol", "_atom_site.Cartn_y", "_atom_site.occupancy"):
    try:
      _cpp_block(_block(_drop_column(atom_site_header + one, tag)))
    except ValueError as e:
      assert tag in str(e), (tag, str(e))
    else:
      raise AssertionError("missing %s accepted" % tag)
  # the label-column assertion keeps its text (iotbx.pdb.utils matches it)
  try:
    _cpp_block(_block(_drop_column(atom_site_header + one, "_atom_site.label_alt_id")))
  except AssertionError as e:
    assert "something is not present" in str(e)
  else:
    raise AssertionError("missing alt_id accepted")
  # F6: no insertion-code column -> icode " " (Python: ""), resid and id_str unchanged
  text = _drop_column(atom_site_header + one, "_atom_site.pdbx_PDB_ins_code")
  ref = _python(_block(text))
  new = _cpp_block(_block(text))
  rg_ref = ref.only_model().only_chain().residue_groups()[0]
  rg_new = new.only_model().only_chain().residue_groups()[0]
  assert rg_ref.icode == "" and rg_new.icode == " "
  assert rg_ref.resid() == rg_new.resid()
  assert ref.atoms()[0].id_str() == new.atoms()[0].id_str()
  # F7: short non-integer seq id -> right-justified (Python: AssertionError);
  # null and blank seq ids rejected (blank ones would merge residue groups)
  text = atom_site_header + "ATOM 1 C C . ALA A 1 1A ? 1.0 2.0 3.0 1.00 10.00 1A ALA A C 1\n"
  try:
    _python(_block(text))
  except AssertionError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; update divergence F7 in iotbx/pdb/mmcif_hierarchy.h")
  h = _cpp_block(_block(text))
  assert h.only_model().only_chain().residue_groups()[0].resseq == "  1A"
  for seq in (".", "?", "''", "' '", "ABCDE"):   # '' and ' ' quoted in the row
    text = atom_site_header + "ATOM 1 C C . ALA A 1 %s ? 1.0 2.0 3.0 1.00 10.00 %s ALA A C 1\n" % (seq, seq)
    try:
      _cpp_block(_block(text))
    except ValueError as e:
      assert "row 1" in str(e), str(e)
      if seq != "ABCDE":
        assert "null value" in str(e), str(e)
    else:
      raise AssertionError("seq id %r accepted" % seq)
  # F9: rule-3 chain break with the same seq id -> new residue group in the new chain
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 0.0 0.0 0.0 1.00 10.00 5 ALA A N 1
ATOM 2 C CA . ALA A 1 5 ? 1.0 0.0 0.0 1.00 10.00 5 ALA A CA 1
ATOM 3 N N . GLY B 2 5 ? 2.0 0.0 0.0 1.00 10.00 5 GLY A N 1
ATOM 4 C CA . GLY B 2 5 ? 3.0 0.0 0.0 1.00 10.00 5 GLY A CA 1
"""
  ref = _python(_block(text))
  assert [len(c.residue_groups()) for c in ref.models()[0].chains()] == [1, 0]
  new = _cpp_block(_block(text))
  assert [len(c.residue_groups()) for c in new.models()[0].chains()] == [1, 1]
  assert [ag.resname for ag in new.models()[0].chains()[0].residue_groups()[0].atom_groups()] == ["ALA"]
  assert [ag.resname for ag in new.models()[0].chains()[1].residue_groups()[0].atom_groups()] == ["GLY"]


def exercise_divergences_f4_f5_f6_output():
  """
  The remaining divergence cases: F6 output identity, F4 without
  _atom_site_anisotrop, F5 on a non-integer atom id.
  """
  one = "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n"
  # F6: without the insertion-code column the icode differs (" " against ""),
  # but make_resid pads it, so show() and the PDB output are identical
  for text in (_drop_column(atom_site_header + one, "_atom_site.pdbx_PDB_ins_code"),
               input_1ab1, input_4edr, input_null_occ_b):
    ref = _python(_block(text))
    new = _cpp_block(_block(text))
    assert [rg.icode for rg in ref.residue_groups()] == [""] * len(list(ref.residue_groups()))
    assert [rg.icode for rg in new.residue_groups()] == [" "] * len(list(new.residue_groups()))
    assert not show_diff(new.as_pdb_string(), ref.as_pdb_string())
    s_ref = StringIO()
    ref.show(out=s_ref)
    s_new = StringIO()
    new.show(out=s_new)
    assert not show_diff(s_new.getvalue(), s_ref.getvalue())
  # F4: no _atom_site.id and no _atom_site_anisotrop -> serial from row order
  two = one + "ATOM 2 C CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1\n"
  text = _drop_column(atom_site_header + two, "_atom_site.id")
  try:
    _python(_block(text))
  except TypeError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; update divergence F4 in iotbx/pdb/mmcif_hierarchy.h")
  h = _cpp_block(_block(text))
  assert [a.serial for a in h.atoms()] == ["    1", "    2"]
  # F5: a non-integer atom id names the tag and the row
  text = atom_site_header + one.replace("ATOM 1 ", "ATOM x1 ")
  try:
    _cpp_block(_block(text))
  except ValueError as e:
    assert "_atom_site.id" in str(e) and "row 1" in str(e), str(e)
  else:
    raise AssertionError("non-integer atom id accepted")


def exercise_single_item_atom_site():
  """
  A one-atom _atom_site given as items rather than a loop. F10: the Python
  builder reads group_PDB, id, pdbx_PDB_model_num, pdbx_PDB_ins_code and
  phenix_scat_dispersion_real/imag without wrapping them into a column, so it
  takes the first character of the value; the new path wraps every column.
  """
  text = "data_t\n"
  tags = [l for l in atom_site_header.splitlines() if l.startswith("_atom_site.")]
  values = "HETATM 15 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1".split()
  for tag, value in zip(tags, values):
    text += "%s %s\n" % (tag, value)
  ref = _python(_block(text))
  if [(a.serial, a.hetero) for a in ref.atoms()] != [("    1", False)]:
    raise AssertionError("the reference builder no longer reads single items by character; update divergence F10 in iotbx/pdb/mmcif_hierarchy.h")
  new = _cpp_block(_block(text))
  assert new.atoms_size() == 1
  assert new.atoms()[0].serial == "   15"
  assert new.atoms()[0].hetero


def exercise_cell_grammar_f11():
  """
  F11: resolution, fp and fdp follow the flex.double grammar and the integer
  cells a strict integer grammar; the Python builder used float() and int().
  """
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.phenix_scat_dispersion_real\n"
    "_atom_site.phenix_scat_dispersion_imag\n_atom_site.phenix_resolution\n")
  row = "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 %s %s %s\n"
  # an esd bracket: 1.5 from C++; float() raised
  text = header + row % ("1.5(3)", "1.5(3)", "1.5(3)")
  try:
    _python(_block(text))
  except ValueError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; update divergence F11 in iotbx/pdb/mmcif_hierarchy.h")
  a = _cpp(text).atoms()[0]
  assert (a.fp, a.fdp, a.resolution) == (1.5, 1.5, 1.5)
  # overflow: ValueError naming the tag from C++; float() gave inf
  tags = ("_atom_site.phenix_scat_dispersion_real",
          "_atom_site.phenix_scat_dispersion_imag", "_atom_site.phenix_resolution")
  for i, tag in enumerate(tags):
    cells = ["0.5", "0.5", "0.5"]
    cells[i] = "1e500"
    text = header + row % tuple(cells)
    a = _python(_block(text)).atoms()[0]
    if not math.isinf((a.fp, a.fdp, a.resolution)[i]):
      raise AssertionError("the reference builder no longer gives inf; update divergence F11 in iotbx/pdb/mmcif_hierarchy.h")
    try:
      _cpp(text)
    except ValueError as e:
      assert tag in str(e) and "row 1" in str(e), str(e)
    else:
      raise AssertionError("overflow in %s accepted" % tag)
  # underscores: int() accepts them (Python 3.6 and later), the C++ grammar does not
  python_int_underscores = sys.version_info >= (3, 6)
  text = atom_site_header + \
    "ATOM 1 C C . ALA A 1 1_0 ? 1.0 2.0 3.0 1.00 10.00 1_0 ALA A C 1\n"
  if python_int_underscores:
    rg = _python(_block(text)).only_model().only_chain().residue_groups()[0]
    if rg.resseq != "  10":
      raise AssertionError("the reference builder no longer reads 1_0 as 10; update divergence F11 in iotbx/pdb/mmcif_hierarchy.h")
  rg = _cpp(text).only_model().only_chain().residue_groups()[0]
  assert rg.resseq == " 1_0", rg.resseq   # not an integer: F7, right-justified
  text = atom_site_header + \
    "ATOM 1_0 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n"
  if python_int_underscores:
    if _python(_block(text)).atoms()[0].serial != "   10":
      raise AssertionError("the reference builder no longer reads 1_0 as 10; update divergence F11 in iotbx/pdb/mmcif_hierarchy.h")
  try:
    _cpp(text)
  except ValueError as e:
    assert "_atom_site.id" in str(e) and "row 1" in str(e), str(e)
  else:
    raise AssertionError("atom id 1_0 accepted")


def exercise_blank_seq_id_f7():
  """F7: a blank seq id is rejected; the Python builder kept a four-blank one."""
  text = atom_site_header + \
    "ATOM 1 C C . ALA A 1 '    ' ? 1.0 2.0 3.0 1.00 10.00 '    ' ALA A C 1\n"
  rg = _python(_block(text)).only_model().only_chain().residue_groups()[0]
  if rg.resseq != "    ":
    raise AssertionError("the reference builder no longer accepts a blank seq id; update divergence F7 in iotbx/pdb/mmcif_hierarchy.h")
  try:
    _cpp(text)
  except ValueError as e:
    assert "row 1" in str(e) and "null value" in str(e), str(e)
  else:
    raise AssertionError("blank seq id accepted")


def exercise_aniso_join_leading_zeros_f1():
  """
  F1: the anisotrop join uses the _atom_site.id string as written; the Python
  builder looked up the stripped hybrid-36 serial, so it joined atom 01 to
  anisotrop 1 and lost the tensor of atom 00001.
  """
  joined = (0.1, 0.2, 0.3, 0.01, 0.02, 0.03)
  unset = (-1.0,) * 6
  for atom_id, aniso_id, python_uij, cpp_uij in (
      ("01", "1", joined, unset),
      ("00001", "00001", unset, joined)):
    text = atom_site_header + \
      "ATOM %s C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n" % atom_id + \
      aniso_u_loop + "%s 0.1 0.2 0.3 0.01 0.02 0.03\n" % aniso_id
    if _python(_block(text)).atoms()[0].uij != python_uij:
      raise AssertionError("the reference builder joins differently now; update divergence F1 in iotbx/pdb/mmcif_hierarchy.h")
    assert _cpp(text).atoms()[0].uij == cpp_uij, (atom_id, aniso_id)


def exercise_reader_uses_cpp():
  """cif_input and iotbx.pdb.input go through the C++ builder."""
  import iotbx.pdb
  text = atom_site_header + "ATOM 1 C C . ALA A 1 1A ? 1.0 2.0 3.0 1.00 10.00 1A ALA A C 1\n"
  h = iotbx.pdb.input(lines=text.splitlines(), source_info=None).construct_hierarchy()
  assert h.only_model().only_chain().residue_groups()[0].resseq == "  1A"   # F7 only the C++ path accepts
  h = mmcif.pdb_hierarchy_builder(_block(text)).hierarchy
  assert h.only_model().only_chain().residue_groups()[0].resseq == "  1A"


def run():
  exercise_cpp_columns_skeleton()
  exercise_cpp_cell_parsing()
  exercise_cpp_core_parity()
  exercise_zero_rows()
  exercise_model_boundaries()
  exercise_chain_break_rule()
  exercise_atom_name_and_element()
  exercise_charge()
  exercise_segid_and_serial()
  exercise_null_rejection_and_nulls_in_b_occ()
  exercise_aniso_join()
  exercise_fp_fdp_resolution_cpp()
  exercise_divergences_f3_f5_f6_f7_f9()
  exercise_divergences_f4_f5_f6_output()
  exercise_single_item_atom_site()
  exercise_cell_grammar_f11()
  exercise_blank_seq_id_f7()
  exercise_aniso_join_leading_zeros_f1()
  exercise_reader_uses_cpp()


if __name__ == "__main__":
  run()
  print("OK")
