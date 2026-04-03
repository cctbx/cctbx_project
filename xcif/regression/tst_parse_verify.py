from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_parse_verify.py
#
# Step 11: Parse verification and edge-case tests.
# Exercises xcif parsing of various CIF constructs (quoted strings,
# semicolon fields, nulls, SU values, scientific notation, etc.)
# and verifies that xcif can re-parse iotbx.cif.model.show() output.
# True round-trip tests (xcif write -> xcif re-parse) are deferred to
# Step 12 when write support is implemented.

import os
import math
import xcif_ext
import iotbx.cif
from iotbx.cif import model as cif_model
import libtbx.load_env
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal

try:
  from six.moves import StringIO
except ImportError:
  from io import StringIO


def exercise_tag_values():
  """Parse -> extract -> verify for simple tag-value pairs."""
  cif_text = (
    "data_rt\n"
    "_cell.length_a 10.500\n"
    "_cell.length_b 20.000\n"
    "_cell.length_c 30.000\n"
    "_cell.angle_alpha 90.00\n"
    "_symmetry.space_group_name_H-M 'P 1'\n"
  )
  doc = xcif_ext.parse(cif_text)
  block = doc[0]

  assert block.name == "rt"
  assert block.find_value("_cell.length_a") == "10.500"
  assert block.find_value("_cell.length_b") == "20.000"
  assert block.find_value("_cell.length_c") == "30.000"
  assert block.find_value("_cell.angle_alpha") == "90.00"
  assert block.find_value("_symmetry.space_group_name_H-M") == "P 1"

  assert approx_equal(xcif_ext.as_double("10.500"), 10.5)
  assert approx_equal(xcif_ext.as_double("90.00"), 90.0)


def exercise_loop_cells():
  """Parse -> extract all columns -> verify each cell."""
  cif_text = (
    "data_loop_rt\n"
    "loop_\n"
    "_atom.id\n"
    "_atom.symbol\n"
    "_atom.x\n"
    "_atom.y\n"
    "_atom.z\n"
    "1 C 1.000 2.000 3.000\n"
    "2 N 4.000 5.000 6.000\n"
    "3 O 7.000 8.000 9.000\n"
  )
  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_atom")
  assert loop is not None

  assert loop.width == 5
  assert loop.length == 3

  # Verify every cell value by row/col
  expected = [
    ["1", "C", "1.000", "2.000", "3.000"],
    ["2", "N", "4.000", "5.000", "6.000"],
    ["3", "O", "7.000", "8.000", "9.000"],
  ]
  for r in range(3):
    for c in range(5):
      actual = loop.value(r, c)
      assert actual == expected[r][c], \
        "Cell (%d,%d): expected '%s', got '%s'" % (r, c, expected[r][c], actual)

  # Verify column extraction
  ids = loop.column_as_flex_int("_atom.id")
  assert list(ids) == [1, 2, 3]

  x = loop.column_as_flex_double("_atom.x")
  assert approx_equal(list(x), [1.0, 4.0, 7.0])

  syms = loop.column_as_flex_string("_atom.symbol")
  assert list(syms) == ["C", "N", "O"]


def exercise_iotbx_show_to_xcif_parse():
  """iotbx.cif model -> show() -> string -> xcif.parse()."""
  # Build a CIF model programmatically via iotbx
  m = cif_model.cif()
  b = cif_model.block()
  b["_cell.length_a"] = "25.000"
  b["_cell.length_b"] = "30.000"
  b["_cell.length_c"] = "35.000"
  b["_symmetry.space_group_name_H-M"] = "P 21 21 21"
  m["roundtrip"] = b

  # Add a loop
  lp = cif_model.loop(header=[
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.Cartn_x",
  ])
  lp.add_row(["1", "C", "10.5"])
  lp.add_row(["2", "N", "20.5"])
  lp.add_row(["3", "O", "30.5"])
  b.add_loop(lp)

  # Serialize to string via iotbx
  sio = StringIO()
  m.show(out=sio)
  cif_string = sio.getvalue()

  # Re-parse with xcif
  doc = xcif_ext.parse(cif_string)
  assert len(doc) == 1, "Expected 1 block, got %d" % len(doc)

  xb = doc[0]
  assert xb.name == "roundtrip"
  assert xb.find_value("_cell.length_a") == "25.000"
  assert xb.find_value("_cell.length_b") == "30.000"
  assert xb.find_value("_symmetry.space_group_name_H-M") == "P 21 21 21"

  xloop = xb.find_loop("_atom_site")
  assert xloop is not None
  assert xloop.width == 3
  assert xloop.length == 3

  x = xloop.column_as_flex_double("_atom_site.Cartn_x")
  assert approx_equal(list(x), [10.5, 20.5, 30.5])

  ids = xloop.column_as_flex_int("_atom_site.id")
  assert list(ids) == [1, 2, 3]


def exercise_null_values():
  """Null/unknown/inapplicable values survive parsing."""
  cif_text = (
    "data_null_rt\n"
    "loop_\n"
    "_val.id\n"
    "_val.x\n"
    "_val.label\n"
    "1 1.5 good\n"
    "2 . ?\n"
    "3 ? .\n"
    "4 3.0 ok\n"
  )
  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_val")

  # Raw string values preserved
  assert loop.value(1, 1) == "."
  assert loop.value(2, 1) == "?"
  assert loop.value(1, 2) == "?"
  assert loop.value(2, 2) == "."

  # Numeric: null -> NaN
  x = loop.column_as_flex_double("_val.x")
  assert approx_equal(x[0], 1.5)
  assert math.isnan(x[1])
  assert math.isnan(x[2])
  assert approx_equal(x[3], 3.0)

  # Predicates
  assert xcif_ext.is_null(".")
  assert xcif_ext.is_null("?")
  assert xcif_ext.is_unknown(".")
  assert not xcif_ext.is_unknown("?")
  assert xcif_ext.is_inapplicable("?")
  assert not xcif_ext.is_inapplicable(".")


def exercise_quoted_strings():
  """Quoted strings have quotes stripped after parsing."""
  cif_text = (
    "data_quote_rt\n"
    "_tag1 'value with spaces'\n"
    "_tag2 \"double quoted\"\n"
    "_tag3 noquotes\n"
    "_tag4 'P 21 21 21'\n"
  )
  doc = xcif_ext.parse(cif_text)
  b = doc[0]

  assert b.find_value("_tag1") == "value with spaces"
  assert b.find_value("_tag2") == "double quoted"
  assert b.find_value("_tag3") == "noquotes"
  assert b.find_value("_tag4") == "P 21 21 21"


def exercise_semicolon_fields():
  """Semicolon text fields survive parsing."""
  cif_text = (
    "data_semi_rt\n"
    "_struct.title\n"
    ";Line one of the title.\n"
    "Line two of the title.\n"
    ";\n"
    "_other.value simple\n"
  )
  doc = xcif_ext.parse(cif_text)
  b = doc[0]

  title = b.find_value("_struct.title")
  assert "Line one" in title, "Missing text: '%s'" % title
  assert "Line two" in title, "Missing text: '%s'" % title
  assert b.find_value("_other.value") == "simple"


def exercise_su_values():
  """Standard uncertainty values survive parsing and extraction."""
  cif_text = (
    "data_su_rt\n"
    "loop_\n"
    "_cell.param\n"
    "_cell.value\n"
    "length_a 50.840(10)\n"
    "length_b 42.770(1)\n"
    "length_c 28.950(2000)\n"
    "angle_alpha 90.00\n"
  )
  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_cell")

  # Raw string values preserved
  vals = loop.column("_cell.value")
  assert vals[0] == "50.840(10)"
  assert vals[1] == "42.770(1)"
  assert vals[2] == "28.950(2000)"
  assert vals[3] == "90.00"

  # as_double strips SU
  assert approx_equal(xcif_ext.as_double("50.840(10)"), 50.840)
  assert approx_equal(xcif_ext.as_double("42.770(1)"), 42.770)
  assert approx_equal(xcif_ext.as_double("28.950(2000)"), 28.950)

  # as_double_with_su extracts both
  v, s = xcif_ext.as_double_with_su("50.840(10)")
  assert approx_equal(v, 50.840)
  assert approx_equal(s, 0.010)

  v, s = xcif_ext.as_double_with_su("42.770(1)")
  assert approx_equal(v, 42.770)
  assert approx_equal(s, 0.001)

  v, s = xcif_ext.as_double_with_su("28.950(2000)")
  assert approx_equal(v, 28.950)
  assert approx_equal(s, 2.000)

  # No SU -> SU is 0
  v, s = xcif_ext.as_double_with_su("90.00")
  assert approx_equal(v, 90.0)
  assert approx_equal(s, 0.0)


def exercise_multi_block():
  """Multi-block documents parsed correctly."""
  cif_text = (
    "data_alpha\n"
    "_tag.a 1\n"
    "data_beta\n"
    "_tag.b 2\n"
    "data_gamma\n"
    "_tag.c 3\n"
  )
  doc = xcif_ext.parse(cif_text)
  assert len(doc) == 3

  assert doc[0].name == "alpha"
  assert doc[1].name == "beta"
  assert doc[2].name == "gamma"

  assert doc[0].find_value("_tag.a") == "1"
  assert doc[1].find_value("_tag.b") == "2"
  assert doc[2].find_value("_tag.c") == "3"

  # find_block case-insensitive
  assert doc.find_block("ALPHA") is not None
  assert doc.find_block("Alpha") is not None
  assert doc.find_block("alpha").find_value("_tag.a") == "1"


def exercise_iotbx_show_with_loop():
  """iotbx build -> show -> xcif parse -> compare with iotbx re-parse."""
  # Build via iotbx
  m = cif_model.cif()
  b = cif_model.block()
  b["_cell.length_a"] = "12.345"
  b["_cell.length_b"] = "23.456"

  lp = cif_model.loop(header=[
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
  ])
  for i in range(10):
    lp.add_row([
      str(i + 1),
      ["C", "N", "O", "S", "P"][i % 5],
      "%.3f" % (i * 1.1),
      "%.3f" % (i * 2.2),
      "%.3f" % (i * 3.3),
    ])
  b.add_loop(lp)
  m["full_rt"] = b

  # Serialize
  sio = StringIO()
  m.show(out=sio)
  cif_string = sio.getvalue()

  # Parse with iotbx (reference) and xcif (test)
  ref_model = iotbx.cif.reader(input_string=cif_string).model()
  doc = xcif_ext.parse(cif_string)

  # Compare
  assert len(doc) == 1
  xb = doc[0]
  rb = list(ref_model.values())[0]

  assert xb.find_value("_cell.length_a") == rb["_cell.length_a"]
  assert xb.find_value("_cell.length_b") == rb["_cell.length_b"]

  xloop = xb.find_loop("_atom_site")
  rloop = rb.get_loop("_atom_site")
  assert xloop.width == rloop.n_columns()
  assert xloop.length == rloop.n_rows()

  # Compare coordinate columns
  for tag in ["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]:
    x_vals = xloop.column_as_flex_double(tag)
    r_vals_str = rloop[tag]
    for j in range(x_vals.size()):
      assert approx_equal(x_vals[j], float(r_vals_str[j])), \
        "%s[%d]: xcif=%f vs ref=%f" % (tag, j, x_vals[j], float(r_vals_str[j]))


def exercise_scientific_notation():
  """Scientific notation values survive parsing."""
  cif_text = (
    "data_sci\n"
    "loop_\n"
    "_val.x\n"
    "1.5e2\n"
    "3.60e-02\n"
    "-1.0E+3\n"
    "0.0\n"
    "1e10\n"
  )
  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_val")
  x = loop.column_as_flex_double("_val.x")

  assert approx_equal(x[0], 150.0)
  assert approx_equal(x[1], 0.036)
  assert approx_equal(x[2], -1000.0)
  assert approx_equal(x[3], 0.0)
  assert approx_equal(x[4], 1e10)


def exercise_empty_block():
  """Blocks with no tags or loops are valid."""
  cif_text = "data_empty\ndata_notempty\n_tag 1\n"
  doc = xcif_ext.parse(cif_text)
  assert len(doc) == 2
  assert doc[0].name == "empty"
  assert doc[1].name == "notempty"
  assert doc[0].find_value("_tag") is None
  assert doc[1].find_value("_tag") == "1"


def exercise_example_cif():
  """Real-world example.cif: parse -> extract -> verify key values."""
  dist_dir = libtbx.env.dist_path("xcif")
  cif_path = os.path.join(dist_dir, "regression", "example.cif")

  doc = xcif_ext.parse_file(cif_path)
  assert len(doc) == 1
  b = doc[0]
  assert b.name.upper() == "1UBQ"

  # Cell parameters with SU
  v, s = xcif_ext.as_double_with_su(b.find_value("_cell.length_a"))
  assert approx_equal(v, 50.840)
  assert approx_equal(s, 0.010)

  v, s = xcif_ext.as_double_with_su(b.find_value("_cell.length_b"))
  assert approx_equal(v, 42.770)
  assert approx_equal(s, 0.001)

  v, s = xcif_ext.as_double_with_su(b.find_value("_cell.length_c"))
  assert approx_equal(v, 28.950)
  assert approx_equal(s, 2.000)

  # Quoted string
  sg = b.find_value("_symmetry.space_group_name_H-M")
  assert sg == "P 21 21 21", "Got '%s'" % sg

  # Semicolon text field
  title = b.find_value("_struct.title")
  assert "ubiquitin" in title.lower(), "Title: '%s'" % title

  # Atom site loop
  loop = b.find_loop("_atom_site")
  assert loop.width == 13
  assert loop.length == 10

  # First atom
  assert loop.value(0, 0) == "1"
  x = loop.column_as_flex_double("_atom_site.Cartn_x")
  assert approx_equal(x[0], 27.340)

  # Reflection stats loop
  rloop = b.find_loop("_reflns")
  assert rloop is not None
  assert rloop.length == 1
  rmerge = rloop.column_as_flex_double("_reflns.pdbx_Rmerge_I_obs")
  assert approx_equal(rmerge[0], 0.036)


if __name__ == "__main__":
  exercise_tag_values()
  exercise_loop_cells()
  exercise_iotbx_show_to_xcif_parse()
  exercise_null_values()
  exercise_quoted_strings()
  exercise_semicolon_fields()
  exercise_su_values()
  exercise_multi_block()
  exercise_iotbx_show_with_loop()
  exercise_scientific_notation()
  exercise_empty_block()
  exercise_example_cif()
  print(format_cpu_times())
  print("OK")
