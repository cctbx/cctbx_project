"""
The C++ mmCIF hierarchy builder (iotbx.pdb.build_mmcif_hierarchy on the columns
of iotbx.pdb.mmcif.atom_site_columns_from_block) against the Python builder
(iotbx.pdb.mmcif._pdb_hierarchy_builder_python): the column binding, cell
parsing, the boundary, naming and field rules, the inputs the Python builder
was fixed for, and the two differences listed in iotbx/pdb/mmcif_hierarchy.h
(error reporting, cell grammar), asserted on both builders.
"""
from six.moves import cStringIO as StringIO
from libtbx.test_utils import show_diff, Exception_expected
import iotbx.cif
from iotbx.cif.builders import CifBuilderError
from iotbx.pdb import mmcif
from iotbx.pdb import atom_site_columns, build_mmcif_hierarchy
from iotbx.pdb import mmcif_parse_cell_double, mmcif_parse_cell_int
from iotbx.pdb.tst_mmcif_hierarchy import input_1ab1, input_4edr, \
  input_null_occ_b, input_single_item_atom_site, atom_site_header, \
  aniso_u_loop, aniso_join_row, aniso_join_uij, aniso_join_cases, \
  seq_id_resseq, seq_id_rejected, drop_column, first_block
from iotbx.regression.tst_mmcif_hierarchy_parity import compare_hierarchies
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
  # a mandatory column of the wrong length is reported by its CIF tag
  cols = atom_site_columns()
  cols.atom_name = flex.std_string(["N", "CA"])
  cols.type_symbol = flex.std_string(["N"])
  try:
    build_mmcif_hierarchy(cols)
  except ValueError as e:
    assert "_atom_site.type_symbol has 1 rows, expected 2" in str(e), str(e)
  else:
    raise Exception_expected
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
    assert "_atom_site.auth_segid has 2 rows, expected 1" in str(e), str(e)
  else:
    raise Exception_expected


def exercise_cpp_cell_parsing():
  # everything flex.double accepts, with the same value bit for bit
  for cell in ("1.5", "-2", "+3.25", ".5", "5.", "1e3", "1E-3", "12.345(6)",
               "1.5+03", "1.5-03", "0", "-0.0", "inf", "nan",
               "1.7976931348623157e308"):
    expected = flex.double(flex.std_string([cell]))[0]
    got = mmcif_parse_cell_double(cell)
    assert repr(got) == repr(expected), (cell, got, expected)
  # cells on which flex.double and this grammar agree on acceptance
  for cell in ("", " 1.5", "1.5 ", "1.5x", "?", ".", "1e999", "1.5(", "(5)",
               "1.5()", "1.5(a)", "abc", "1,5", "1.5d0", "1.5(99999999999)",
               "1.5(2147483647)", "1.5(2147483648)", "1.5(+5)", "1.5(-5)",
               "1.5\x00abc", "1e", "1E", "1+", "1-", "\xa01.5", "1.5\xa0"):
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
  # rejected here whatever the platform's flex.double does: hex floats
  # (libc++ reads them), underflow and subnormal values (libstdc++ reads the
  # subnormals); the message quotes the cell
  for cell in ("0x10", "0X10", "0x1e", "0x1p-2", "1e-400", "-1e-400",
               "4.9e-324", "1e-310", "1e-310(2)"):
    try:
      mmcif_parse_cell_double(cell)
    except ValueError as e:
      assert '"%s"' % cell in str(e), str(e)
      assert ("value part" in str(e)) == ("(" in cell), str(e)
    else:
      raise AssertionError("accepted %r" % cell)
  assert mmcif_parse_cell_int("12") == 12
  assert mmcif_parse_cell_int(" -7 ") == -7
  assert mmcif_parse_cell_int("+5") == 5
  assert mmcif_parse_cell_int("\t12\n") == 12
  # whitespace is the ASCII set whatever the locale: a no-break space is not
  # trimmed (as a Python str, UTF-8 encoded; as bytes, the Latin-1 byte)
  for cell in ("", "1.5", "1a", "a1", "?", ".", "1 2", "--1",
               "99999999999999999999", "12\xa0", "\xa012",
               b"12\xa0", b"\xa012"):
    try:
      mmcif_parse_cell_int(cell)
    except ValueError as e:
      assert str(e).isascii(), repr(str(e))
    else:
      raise AssertionError("accepted %r" % cell)


# column order of atom_site_header rows:
# group id elem name alt comp lasym ent lseq icode x y z occ b
# aseq acomp aasym aname model


def _python(block):
  """The Python builder, the reference."""
  return mmcif._pdb_hierarchy_builder_python(block).hierarchy


def _cpp_block(block):
  """The C++ builder called directly on the columns of
  atom_site_columns_from_block."""
  return build_mmcif_hierarchy(mmcif.atom_site_columns_from_block(block))


def _cpp(text):
  return _cpp_block(first_block(text))


def _both(text):
  """(Python build, C++ build), both unsorted."""
  return (_python(first_block(text)), _cpp(text))


def _same(text):
  """Both builders build text and the hierarchies compare equal; returns the
  C++ one."""
  ref, new = _both(text)
  compare_hierarchies(ref, new)
  return new


def _both_raise(text, fragments, python_error=Exception):
  """Both builders reject text: the Python one with python_error, the C++
  one with a ValueError whose message is ASCII and contains every
  fragment."""
  block = first_block(text)
  try:
    _python(block)
  except python_error:
    pass
  else:
    raise Exception_expected
  try:
    _cpp_block(block)
  except ValueError as e:
    assert str(e).isascii(), repr(str(e))
    for fragment in fragments:
      assert fragment in str(e), (fragment, str(e))
  else:
    raise Exception_expected


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
  new = _same(text)
  assert new.models_size() == 2
  assert [c.id for c in new.models()[0].chains()] == ["A", "A", "B"]
  rgs = new.models()[0].chains()[0].residue_groups()
  assert [rg.resid() for rg in rgs] == ["   5 ", "   6 ", "   6A"]
  assert [ag.altloc for ag in rgs[0].atom_groups()] == ["", "A", "B"]
  assert new.models()[0].chains()[0].atoms()[0].name == " N  "
  assert new.models()[0].chains()[1].atoms()[1].name == "FE  "
  assert new.models()[0].chains()[1].atoms()[0].hetero
  # two blank-altloc residue names in one residue group: each new blank group
  # is inserted at the front, so the order is reversed (as in the Python
  # builder)
  text = atom_site_header + """\
HETATM 1 O O . HOH A 1 . ? 1.0 2.0 3.0 1.00 10.00 7 HOH A O 1
HETATM 2 NA NA . NA A 1 . ? 1.0 2.0 3.0 1.00 10.00 7 NA A NA 1
"""
  new = _same(text)
  rg0 = list(new.residue_groups())[0]   # root.residue_groups() is a generator
  assert [ag.resname for ag in rg0.atom_groups()] == [" NA", "HOH"]


def exercise_zero_rows():
  # _atom_site present with zero rows -> empty hierarchy, also when
  # _atom_site_anisotrop has rows (there is no atom to join them to)
  for text in (atom_site_header,
               atom_site_header + aniso_u_loop + aniso_join_row % "1"):
    new = _same(text)
    assert new.models_size() == 0


def exercise_model_boundaries():
  # a new model starts wherever the model number changes, so non-adjacent
  # repeats of a number are separate models
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 1
ATOM 2 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 2
ATOM 3 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 1
"""
  new = _same(text)
  assert [m.id for m in new.models()] == ["1", "2", "1"]
  # single model: id becomes ""
  new = _same(atom_site_header +
    "ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 7\n")
  assert new.models()[0].id == ""
  # no model number column: one model with id ""
  new = _same(drop_column(atom_site_header +
    "ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 5 ALA A N 7\n",
    "_atom_site.pdbx_PDB_model_num"))
  assert new.models_size() == 1 and new.models()[0].id == ""


def exercise_chain_break_rule():
  """
  The chain rules: a new chain where the model or auth_asym_id changes, or
  where label_asym_id changes after a polymer residue (rule 3); a new chain
  always starts a residue group.
  """
  # polymer, then ligand and water in the same auth chain: split after the
  # polymer only
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
ATOM 2 N N . GLY A 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1
HETATM 3 C C1 . LIG B 2 . ? 1.0 2.0 3.0 1.00 10.00 101 LIG A C1 1
HETATM 4 C C1 . LIG C 2 . ? 1.0 2.0 3.0 1.00 10.00 102 LIG A C1 1
HETATM 5 O O . HOH D 3 . ? 1.0 2.0 3.0 1.00 10.00 201 HOH A O 1
HETATM 6 O O . HOH D 3 . ? 1.0 2.0 3.0 1.00 10.00 202 HOH A O 1
"""
  new = _same(text)
  assert [len(c.residue_groups()) for c in new.models()[0].chains()] == [2, 4]
  # two same-id protein chains as the Phenix writer emits them: split
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
ATOM 2 N N . GLY A 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1
ATOM 3 N N . ALA B 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1
ATOM 4 N N . GLY B 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1
"""
  new = _same(text)
  assert [c.id for c in new.models()[0].chains()] == ["A", "A"]
  # a D-amino acid before the label change: no split today, reproduced
  text = atom_site_header + """\
ATOM 1 N N . DAL A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 DAL A N 1
ATOM 2 N N . ALA B 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 ALA A N 1
"""
  new = _same(text)
  assert len(new.models()[0].chains()) == 1
  # a label_asym_id break with the sequence id and insertion code unchanged:
  # the new chain gets its own residue group
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 0.0 0.0 0.0 1.00 10.00 5 ALA A N 1
ATOM 2 C CA . ALA A 1 5 ? 1.0 0.0 0.0 1.00 10.00 5 ALA A CA 1
ATOM 3 N N . GLY B 2 5 ? 2.0 0.0 0.0 1.00 10.00 5 GLY A N 1
ATOM 4 C CA . GLY B 2 5 ? 3.0 0.0 0.0 1.00 10.00 5 GLY A CA 1
"""
  new = _same(text)
  chains = new.models()[0].chains()
  assert [len(c.residue_groups()) for c in chains] == [1, 1]
  assert [[ag.resname for ag in c.residue_groups()[0].atom_groups()]
          for c in chains] == [["ALA"], ["GLY"]]
  # auth_break clears link_to_previous on a later residue group only
  text = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.auth_break\n") + """\
ATOM 1 N N . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A N 1 1
ATOM 2 N N . GLY A 1 2 ? 1.0 2.0 3.0 1.00 10.00 2 GLY A N 1 0
ATOM 3 N N . SER A 1 3 ? 1.0 2.0 3.0 1.00 10.00 3 SER A N 1 1
"""
  new = _same(text)
  assert [rg.link_to_previous for rg in new.residue_groups()] == [
    True, True, False]


def exercise_null_chain_id():
  """A null auth_asym_id: an assertion in the Python builder, an error naming
  tag and row from the C++ builder (exercise_reader_uses_cpp has it as a
  CifBuilderError through the reader)."""
  one = "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA %s C 1\n"
  for chain_id in (".", "?", "' '"):
    text = atom_site_header + one % chain_id
    try:
      _python(first_block(text))
    except AssertionError as e:
      assert "empty auth_asym_id" in str(e), str(e)
    else:
      raise AssertionError(
        "chain id %r accepted by the Python builder" % chain_id)
    try:
      _cpp(text)
    except ValueError as e:
      assert "_atom_site.auth_asym_id" in str(e) and "row 1" in str(e), str(e)
      assert "null" in str(e), str(e)
    else:
      raise AssertionError(
        "chain id %r accepted by the C++ builder" % chain_id)


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
  new = _same(text)
  assert [a.name for a in new.atoms()] == [
    " CA ", "FE1 ", "HD21", " C  ", "CA  "]
  # over-long fields: both builders raise, the C++ one naming tag and row
  for bad_row, tag in (
      ("ATOM 1 C HD211 . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A HD211 1",
       "label_atom_id"),
      ("ATOM 1 FEX FE . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A FE 1",
       "_atom_site.type_symbol"),
      ("ATOM 1 C C . ALA A 1 1 AB 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1",
       "pdbx_PDB_ins_code")):
    _both_raise(atom_site_header + bad_row + "\n", (tag, "row 1"))
  # a multi-character altloc is accepted by both builders (atom_group altloc
  # is a std::string, not a fixed-width field)
  new = _same(atom_site_header +
    "ATOM 1 C C AB ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n")
  assert list(new.residue_groups())[0].atom_groups()[0].altloc == "AB"
  # a blank name with a one- or two-letter element: an error from both
  # builders (IndexError in format_pdb_atom_name); with a blank element it
  # is four blanks from both, which is what the writer emits for a PDB atom
  # with blank name and element columns (exercise_writer_round_trip)
  row = "ATOM 1 %s %s . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A ' ' 1\n"
  for element in ("C", "FE"):
    _both_raise(atom_site_header + row % (element, "' '"),
      ("label_atom_id", "row 1", "blank"), python_error=IndexError)
  new = _same(atom_site_header + row % ("''", "' '"))
  assert new.atoms()[0].name == "    " and new.atoms()[0].element == ""


def exercise_non_ascii_cells():
  """
  A text cell holds CIF 1.1 text: printable ASCII, tab, line feed and
  carriage return. The C++ builder rejects any other byte in every text
  column it reads, with the byte's position, in an ASCII message. The
  Python builder raised UnicodeDecodeError on a byte it could not decode
  and built valid UTF-8 text (iotbx/pdb/mmcif_hierarchy.h).
  """
  row = "ATOM 1 C %s . %s A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 %s %s C 1\n"
  # a tab inside a quoted cell and the line feed of a semicolon text field
  # are CIF text: both builders build them, the name stripped as always
  new = _same(atom_site_header + row % ("'A\tB'", "ALA", "ALA", "A"))
  assert new.atoms()[0].name == "A\tB "
  new = _same(atom_site_header + row % ("\n;CA\n;\n", "ALA", "ALA", "A"))
  assert new.atoms()[0].name == " CA "
  new = _same(atom_site_header + row % ("CA", "ALA", "ALA", "\n;A\n;\n"))
  assert [c.id for c in new.only_model().chains()] == ["A\n"]
  # as a Python str a non-ASCII cell reaches C++ as UTF-8, two bytes here:
  # the Python builder builds the resname (a std::string), the C++ builder
  # rejects it
  text = atom_site_header + row % ("C", "AL\xe9", "AL\xe9", "A")
  h = _python(first_block(text))
  assert list(h.residue_groups())[0].atom_groups()[0].resname == "AL\xe9"
  try:
    _cpp(text)
  except ValueError as e:
    assert str(e).isascii(), repr(str(e))
    assert "comp_id" in str(e) and "row 1" in str(e), str(e)
    assert "byte 3" in str(e), str(e)
    assert '"AL\\xc3\\xa9"' in str(e), str(e)
  else:
    raise Exception_expected
  # a segid cut to four code points in Python overflows the four-byte field:
  # both raise
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.auth_segid\n")
  _both_raise(header +
    "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 ABC\xe9\n",
    ("auth_segid", "row 1", "byte 4"), python_error=ValueError)
  # every text column read is checked, label_asym_id (compared, never
  # stored) included, before any length check; the integer grammar of the
  # other columns quotes the cell the same way. A lone Latin-1 byte, which
  # the Python builder cannot decode, and a NUL.
  cols = mmcif.atom_site_columns_from_block(
    first_block(atom_site_header + row % ("C", "ALA", "ALA", "A")))
  cif = "outside the CIF 1.1 character set"
  for member, tag, reason in (
      ("type_symbol", "type_symbol", cif), ("atom_name", "atom_id", cif),
      ("alt_id", "alt_id", cif), ("label_asym_id", "label_asym_id", cif),
      ("auth_asym_id", "auth_asym_id", cif), ("comp_id", "comp_id", cif),
      ("model_num", "model_num", cif), ("ins_code", "ins_code", cif),
      ("seq_id", "auth_seq_id", "hybrid-36"),
      ("atom_id", "_atom_site.id", "hybrid-36 range"),
      ("formal_charge", "pdbx_formal_charge", "not an integer charge")):
    good = getattr(cols, member)
    for cell, quoted in ((b"\xe9", '"\\xe9"'), (b"A\x00", '"A\\x00"')):
      setattr(cols, member, flex.std_string([cell]))
      try:
        build_mmcif_hierarchy(cols)
      except ValueError as e:
        assert str(e).isascii(), repr(str(e))
        assert tag in str(e) and "row 1" in str(e), (member, str(e))
        assert reason in str(e) and quoted in str(e), (member, str(e))
      else:
        raise AssertionError("%r in %s accepted" % (cell, member))
    setattr(cols, member, good)
  build_mmcif_hierarchy(cols)   # restored: builds


def exercise_column_fallbacks():
  """
  The column fallbacks of atom_site_columns_from_block, through both
  builders: auth_atom_id for label_atom_id, label_comp_id for auth_comp_id,
  label_seq_id for auth_seq_id, and each asym id for the other.
  """
  # label: ALA A 5 N/CA; auth: GLY B 7 NX/CAX
  text = atom_site_header + """\
ATOM 1 N N . ALA A 1 5 ? 1.0 2.0 3.0 1.00 10.00 7 GLY B NX 1
ATOM 2 C CA . ALA A 1 5 ? 1.5 2.0 3.0 1.00 11.00 7 GLY B CAX 1
"""
  new = _same(text)
  assert [a.name for a in new.atoms()] == [" N  ", " CA "]
  assert [c.id for c in new.only_model().chains()] == ["B"]
  new = _same(drop_column(text, "_atom_site.label_atom_id"))
  assert [a.name for a in new.atoms()] == [" NX ", " CAX"]
  new = _same(drop_column(text, "_atom_site.auth_comp_id"))
  assert list(new.residue_groups())[0].atom_groups()[0].resname == "ALA"
  new = _same(drop_column(text, "_atom_site.auth_seq_id"))
  assert [rg.resseq for rg in new.residue_groups()] == ["   5"]
  new = _same(drop_column(text, "_atom_site.auth_asym_id"))
  assert [c.id for c in new.only_model().chains()] == ["A"]
  new = _same(drop_column(text, "_atom_site.label_asym_id"))
  assert [c.id for c in new.only_model().chains()] == ["B"]


def exercise_charge():
  # formal charge: the sign comes from a leading or trailing '-', zero has no
  # sign, ? and . leave the charge blank
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.pdbx_formal_charge\n")
  row = "ATOM %s C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 %s\n"
  rows = ""
  for i, charge in enumerate(("2+", "+2", "-1", "1-", "0", "?", ".")):
    rows += row % (i + 1, charge)
  new = _same(header + rows)
  assert [a.charge for a in new.atoms()] == [
    "2+", "2+", "1-", "1-", "0", "", ""]
  row = row % (1, "%s")
  for charge in ("x+", "10+"):
    _both_raise(header + row % charge, ("pdbx_formal_charge", "row 1"),
      python_error=ValueError)


def exercise_segid_and_serial():
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.auth_segid\n")
  text = header + """\
ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 AB
ATOM 100000 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 ABCDEF
"""
  new = _same(text)
  assert [a.segid for a in new.atoms()] == ["AB  ", "ABCD"]
  assert [a.serial for a in new.atoms()] == ["    1", "A0000"]
  # no _atom_site.id: serials from the row order
  text = drop_column(atom_site_header + """\
ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1
ATOM 2 C CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1
""", "_atom_site.id")
  new = _same(text)
  assert [a.serial for a in new.atoms()] == ["    1", "    2"]
  # a non-integer atom id: both raise, the C++ one naming the tag and the row
  _both_raise(atom_site_header +
    "ATOM x1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n",
    ("_atom_site.id", "row 1"), python_error=ValueError)


def exercise_seq_ids():
  """
  resseq: an integer encoded into the four characters (hybrid-36 above 9999),
  or a four-character hybrid-36 resseq as cctbx wrote until June 2023;
  anything else, and an integer that does not fit, is an error from both
  builders (CifBuilderError in Python, located ValueError from C++).
  """
  row = "ATOM 1 C C . ALA A 1 %s ? 1.0 2.0 3.0 1.00 10.00 %s ALA A C 1\n"
  for seq_id, resseq in seq_id_resseq:
    new = _same(atom_site_header + row % (seq_id, seq_id))
    rg = new.only_model().only_chain().residue_groups()[0]
    assert rg.resseq == resseq, (seq_id, rg.resseq)
  for seq_id, reason in seq_id_rejected:
    _both_raise(atom_site_header + row % (seq_id, seq_id),
      ("_atom_site.auth_seq_id", "row 1", reason),
      python_error=CifBuilderError)


def exercise_missing_columns():
  one = "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n"
  # no insertion-code column: icode " " as with a "?" column; show() and the
  # PDB output of the fixtures without the column are the same from both
  new = _same(drop_column(atom_site_header + one,
    "_atom_site.pdbx_PDB_ins_code"))
  assert [rg.icode for rg in new.residue_groups()] == [" "]
  for text in (input_1ab1, input_4edr, input_null_occ_b):
    ref, new = _both(text)
    compare_hierarchies(ref, new)
    icodes = [rg.icode for rg in new.residue_groups()]
    assert icodes == [" "] * len(icodes)
    assert not show_diff(new.as_pdb_string(), ref.as_pdb_string())
    s_ref = StringIO()
    ref.show(out=s_ref)
    s_new = StringIO()
    new.show(out=s_new)
    assert not show_diff(s_new.getvalue(), s_ref.getvalue())
  # a missing type_symbol or numeric column is a CifBuilderError naming
  # the tag from atom_site_columns_from_block (the Python builder raises too,
  # without the tag)
  for tag in ("_atom_site.type_symbol", "_atom_site.Cartn_y",
              "_atom_site.occupancy"):
    text = drop_column(atom_site_header + one, tag)
    try:
      _python(first_block(text))
    except Exception:
      pass
    else:
      raise AssertionError("missing %s accepted by the Python builder" % tag)
    try:
      _cpp_block(first_block(text))
    except CifBuilderError as e:
      assert tag in str(e), (tag, str(e))
    else:
      raise AssertionError("missing %s accepted" % tag)
  # the label-column assertion keeps its text (iotbx.pdb.utils matches it)
  text = drop_column(atom_site_header + one, "_atom_site.label_alt_id")
  for build in (_python, _cpp_block):
    try:
      build(first_block(text))
    except AssertionError as e:
      assert "something is not present" in str(e)
    else:
      raise Exception_expected


def exercise_null_rejection_and_nulls_in_b_occ():
  # occupancy and B accept . and ? as 0.0 (integrative models write them);
  # coordinates do not, and a malformed number is rejected with its location
  text = atom_site_header + \
    "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 ? . 1 ALA A C 1\n"
  new = _same(text)
  assert new.atoms()[0].occ == 0.0 and new.atoms()[0].b == 0.0
  row = "ATOM 1 C C . ALA A 1 1 ? %s 2.0 3.0 1.0 %s 1 ALA A C 1\n"
  _both_raise(atom_site_header + row % ("?", "1.0"),
    ("_atom_site.Cartn_x row 1", 'Invalid floating-point value: "?"'),
    python_error=ValueError)
  _both_raise(atom_site_header + row % ("1.0", "1.0.1"),
    ("_atom_site.B_iso_or_equiv row 1",
     'Invalid floating-point value: "1.0.1"'), python_error=ValueError)
  # the cell is quoted once, inside the grammar's own message
  try:
    _cpp(atom_site_header + row % ("x", "1.0"))
  except ValueError as e:
    assert str(e) == \
      '_atom_site.Cartn_x row 1: Invalid floating-point value: "x"', str(e)
  else:
    raise Exception_expected


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
  new = _same(text)
  u = [a.uij for a in new.atoms()]
  assert u[0] == (0.1, 0.2, 0.3, 0.01, 0.02, 0.03)
  assert u[1] == (0.7, 0.8, 0.9, 0.07, 0.08, 0.09)
  assert u[2] == (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0)
  # the join is on the _atom_site.id string as written, then on the hybrid-36
  # serial (the cases of tst_mmcif_hierarchy.py)
  for atom_id, aniso_id, joined in aniso_join_cases:
    text = atom_site_header + \
      "ATOM %s C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n" \
      % atom_id + aniso_u_loop + aniso_join_row % aniso_id
    new = _same(text)
    expected = aniso_join_uij if joined else (-1.0,) * 6
    assert new.atoms()[0].uij == expected, (atom_id, aniso_id)
  # B tensors are converted to U, bit for bit the same in both builders
  textb = atom_site_header + rows[:rows.index("ATOM 2")] + \
    aniso_u_loop.replace("U[", "B[") + "1 8.0 8.0 8.0 0.0 0.0 0.0\n"
  new = _same(textb)
  expected = adptbx.b_as_u((8.0, 8.0, 8.0, 0.0, 0.0, 0.0))
  assert new.atoms()[0].uij == expected, (new.atoms()[0].uij, expected)
  # no _atom_site.id with anisotrop rows to join: an error from both
  text4 = drop_column(atom_site_header + rows[:rows.index("ATOM 2")],
    "_atom_site.id") + aniso_u_loop + aniso_join_row % "1"
  _both_raise(text4, ("_atom_site.id",), python_error=CifBuilderError)
  # a bad aniso cell names the anisotrop row
  textbad = atom_site_header + rows + aniso_u_loop + \
    "1 0.1 0.2 x 0.01 0.02 0.03\n"
  _both_raise(textbad, ("_atom_site_anisotrop.U[3][3]", "row 1", '"x"'),
    python_error=CifBuilderError)


def exercise_fp_fdp_resolution_cpp():
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.phenix_scat_dispersion_real\n"
    "_atom_site.phenix_scat_dispersion_imag\n_atom_site.phenix_resolution\n")
  text = header + """\
ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 0.5 1.5 2.0
ATOM 2 C CA . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A CA 1 ? . ?
"""
  new = _same(text)
  a = new.atoms()
  assert (a[0].fp, a[0].fdp, a[0].resolution) == (0.5, 1.5, 2.0)
  assert (a[1].fp, a[1].fdp, a[1].resolution) == (0.0, 0.0, 0.0)


def exercise_single_item_atom_site():
  """A one-atom _atom_site given as items rather than a loop: every item is
  read whole by both builders."""
  new = _same(input_single_item_atom_site)
  assert new.atoms_size() == 1
  a = new.atoms()[0]
  assert a.serial == "   15" and a.hetero
  assert (a.fp, a.fdp) == (0.5, 1.5)
  assert a.parent().parent().icode == "B"


def exercise_cell_grammar():
  """
  resolution, fp and fdp follow the flex.double grammar and the integer
  cells a strict integer grammar; the Python builder uses float() and int().
  """
  header = atom_site_header.replace("_atom_site.pdbx_PDB_model_num\n",
    "_atom_site.pdbx_PDB_model_num\n_atom_site.phenix_scat_dispersion_real\n"
    "_atom_site.phenix_scat_dispersion_imag\n_atom_site.phenix_resolution\n")
  row = "ATOM 1 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1 %s %s %s\n"
  # an esd bracket: 1.5 from C++; float() raises
  text = header + row % ("1.5(3)", "1.5(3)", "1.5(3)")
  try:
    _python(first_block(text))
  except ValueError:
    pass
  else:
    raise AssertionError("the reference builder no longer raises; "
      "update the cell-grammar note in iotbx/pdb/mmcif_hierarchy.h")
  a = _cpp(text).atoms()[0]
  assert (a.fp, a.fdp, a.resolution) == (1.5, 1.5, 1.5)
  # an exponent without its E: 1500 from C++; float() raises
  text = header + row % ("1.5+03", "1.5-03", "1.5+03")
  try:
    _python(first_block(text))
  except ValueError:
    pass
  else:
    raise AssertionError("the reference builder reads a bare exponent now; "
      "update the cell-grammar note in iotbx/pdb/mmcif_hierarchy.h")
  a = _cpp(text).atoms()[0]
  assert (a.fp, a.fdp, a.resolution) == (1500.0, 0.0015, 1500.0)
  # blanks around the number: float() strips them, the grammar does not
  text = header + row % ("' 1.5'", "'1.5 '", "1.5")
  a = _python(first_block(text)).atoms()[0]
  assert (a.fp, a.fdp) == (1.5, 1.5)
  try:
    _cpp(text)
  except ValueError as e:
    assert "phenix_scat_dispersion_real" in str(e), str(e)
    assert '" 1.5"' in str(e), str(e)
  else:
    raise Exception_expected
  # overflow, underflow and a subnormal: ValueError naming the tag from C++;
  # float() gave inf, 0.0 and the subnormal
  tags = ("_atom_site.phenix_scat_dispersion_real",
          "_atom_site.phenix_scat_dispersion_imag",
          "_atom_site.phenix_resolution")
  for cell, python_value in (("1e500", float("inf")), ("1e-400", 0.0),
                             ("4.9e-324", 5e-324)):
    for i, tag in enumerate(tags):
      cells = ["0.5", "0.5", "0.5"]
      cells[i] = cell
      text = header + row % tuple(cells)
      a = _python(first_block(text)).atoms()[0]
      if (a.fp, a.fdp, a.resolution)[i] != python_value:
        raise AssertionError("the reference builder reads %s differently "
          "now; update the cell-grammar note in iotbx/pdb/mmcif_hierarchy.h"
          % cell)
      try:
        _cpp(text)
      except ValueError as e:
        assert tag in str(e) and "row 1" in str(e), str(e)
        assert '"%s"' % cell in str(e), str(e)
      else:
        raise AssertionError("%s in %s accepted" % (cell, tag))
  # underscores: int() accepts them, the C++ grammar does not, and 1_0 is
  # no four-character hybrid-36 resseq either
  text = atom_site_header + \
    "ATOM 1 C C . ALA A 1 1_0 ? 1.0 2.0 3.0 1.00 10.00 1_0 ALA A C 1\n"
  rg = _python(first_block(text)).only_model().only_chain().residue_groups()[0]
  if rg.resseq != "  10":
    raise AssertionError("the reference builder no longer reads 1_0 as "
      "10; update the cell-grammar note in iotbx/pdb/mmcif_hierarchy.h")
  try:
    _cpp(text)
  except ValueError as e:
    assert "_atom_site.auth_seq_id" in str(e) and "row 1" in str(e), str(e)
  else:
    raise Exception_expected
  text = atom_site_header + \
    "ATOM 1_0 C C . ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1 ALA A C 1\n"
  if _python(first_block(text)).atoms()[0].serial != "   10":
    raise AssertionError("the reference builder no longer reads 1_0 as "
      "10; update the cell-grammar note in iotbx/pdb/mmcif_hierarchy.h")
  try:
    _cpp(text)
  except ValueError as e:
    assert "_atom_site.id" in str(e) and "row 1" in str(e), str(e)
  else:
    raise Exception_expected


def run():
  exercise_cpp_columns_skeleton()
  exercise_cpp_cell_parsing()
  exercise_cpp_core_parity()
  exercise_zero_rows()
  exercise_model_boundaries()
  exercise_chain_break_rule()
  exercise_null_chain_id()
  exercise_atom_name_and_element()
  exercise_non_ascii_cells()
  exercise_column_fallbacks()
  exercise_charge()
  exercise_segid_and_serial()
  exercise_seq_ids()
  exercise_missing_columns()
  exercise_null_rejection_and_nulls_in_b_occ()
  exercise_aniso_join()
  exercise_fp_fdp_resolution_cpp()
  exercise_single_item_atom_site()
  exercise_cell_grammar()


if __name__ == "__main__":
  run()
  print("OK")
