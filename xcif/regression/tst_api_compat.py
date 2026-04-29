from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_api_compat.py
#
# Step 11: API compatibility tests.
# Parses identical CIF inputs through both xcif_ext and iotbx.cif.reader,
# then verifies that the extracted data matches.

import os
import math
import xcif
import xcif_ext
import iotbx.cif
import libtbx.load_env
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal

# ---------------------------------------------------------------------------
# Shared CIF strings used by multiple tests
# ---------------------------------------------------------------------------

SIMPLE_CIF = """\
data_compat
_cell.length_a 10.500
_cell.length_b 20.000
_cell.length_c 30.000
_cell.angle_alpha 90.00
_cell.angle_beta 90.00
_cell.angle_gamma 90.00
_symmetry.space_group_name_H-M 'P 21 21 21'
"""

LOOP_CIF = """\
data_loop_test
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_comp_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
1 N MET 27.340 24.430 2.614 1.00 9.67
2 C MET 26.266 25.413 2.842 1.00 10.38
3 C MET 26.913 26.639 3.531 1.00 9.40
4 O MET 27.886 26.463 4.263 1.00 8.32
5 C GLN 25.112 24.880 3.649 1.00 14.45
"""

MULTI_BLOCK_CIF = """\
data_block_A
_cell.length_a 10.0
_entry.id A

data_block_B
_cell.length_a 20.0
_entry.id B

data_block_C
_cell.length_a 30.0
_entry.id C
"""

NULLS_CIF = """\
data_nulls
loop_
_val.id
_val.x
_val.note
1 1.500 good
2 . ?
3 ? .
4 3.000 ok
"""

QUOTED_CIF = """\
data_quoted
_tag1 'single quoted value'
_tag2 "double quoted value"
_tag3 unquoted
"""

SEMICOLON_CIF = """\
data_semi
_struct.title
;This is a multi-line
text field value.
;
_other.tag simple
"""

SAVE_FRAME_CIF = """\
data_dict
_dict.title test_dictionary
save_my_category
_category.id my_category
_category.mandatory_code no
save_
"""


def exercise_block_names():
  """Both parsers find the same block names (case-preserved)."""
  doc = xcif_ext.parse(MULTI_BLOCK_CIF)
  cif_model = iotbx.cif.reader(input_string=MULTI_BLOCK_CIF).model()

  xcif_names = [doc[i].name for i in range(len(doc))]
  iotbx_names = list(cif_model.keys())
  assert xcif_names == iotbx_names, \
    "Block names differ: xcif=%s vs iotbx=%s" % (xcif_names, iotbx_names)


def exercise_tag_values():
  """Single tag-value pairs match between parsers."""
  doc = xcif_ext.parse(SIMPLE_CIF)
  cif_model = iotbx.cif.reader(input_string=SIMPLE_CIF).model()

  xb = doc[0]
  ib = list(cif_model.values())[0]

  tags_to_check = [
    "_cell.length_a",
    "_cell.length_b",
    "_cell.length_c",
    "_cell.angle_alpha",
    "_symmetry.space_group_name_H-M",
  ]
  for tag in tags_to_check:
    xval = xb.find_value(tag)
    ival = ib[tag]
    assert xval == ival, \
      "Tag %s: xcif='%s' vs iotbx='%s'" % (tag, xval, ival)


def exercise_tag_case_insensitive():
  """Case-insensitive tag lookup works in both parsers."""
  doc = xcif_ext.parse(SIMPLE_CIF)
  cif_model = iotbx.cif.reader(input_string=SIMPLE_CIF).model()

  xb = doc[0]
  ib = list(cif_model.values())[0]

  # xcif
  assert xb.find_value("_CELL.LENGTH_A") == "10.500"
  assert xb.find_value("_cell.LENGTH_A") == "10.500"
  # iotbx
  assert ib["_CELL.LENGTH_A"] == "10.500"
  assert ib["_cell.LENGTH_A"] == "10.500"


def exercise_loop_columns():
  """Loop column extraction matches between parsers."""
  doc = xcif_ext.parse(LOOP_CIF)
  cif_model = iotbx.cif.reader(input_string=LOOP_CIF).model()

  xloop = doc[0].find_loop("_atom_site")
  iblock = list(cif_model.values())[0]

  # Compare string columns
  x_syms = xloop.column("_atom_site.type_symbol")
  i_syms = list(iblock["_atom_site.type_symbol"])
  assert x_syms == i_syms, \
    "type_symbol: xcif=%s vs iotbx=%s" % (x_syms, i_syms)

  x_comp = xloop.column("_atom_site.label_comp_id")
  i_comp = list(iblock["_atom_site.label_comp_id"])
  assert x_comp == i_comp, \
    "label_comp_id: xcif=%s vs iotbx=%s" % (x_comp, i_comp)

  # Compare numeric columns via flex
  x_cartn = xloop.column_as_flex_double("_atom_site.Cartn_x")
  i_cartn_str = iblock["_atom_site.Cartn_x"]
  for j in range(x_cartn.size()):
    ival = float(i_cartn_str[j])
    assert approx_equal(x_cartn[j], ival), \
      "Cartn_x[%d]: xcif=%f vs iotbx=%f" % (j, x_cartn[j], ival)

  # Compare integer column
  x_ids = xloop.column_as_flex_int("_atom_site.id")
  i_ids_str = iblock["_atom_site.id"]
  for j in range(x_ids.size()):
    assert x_ids[j] == int(i_ids_str[j]), \
      "id[%d]: xcif=%d vs iotbx=%d" % (j, x_ids[j], int(i_ids_str[j]))


def exercise_loop_dimensions():
  """Loop width and length match between parsers."""
  doc = xcif_ext.parse(LOOP_CIF)
  cif_model = iotbx.cif.reader(input_string=LOOP_CIF).model()

  xloop = doc[0].find_loop("_atom_site")
  iloop = list(cif_model.values())[0].get_loop("_atom_site")

  assert xloop.width == iloop.n_columns(), \
    "width: xcif=%d vs iotbx=%d" % (xloop.width, iloop.n_columns())
  assert xloop.length == iloop.n_rows(), \
    "length: xcif=%d vs iotbx=%d" % (xloop.length, iloop.n_rows())


def exercise_null_handling():
  """Null/unknown values handled identically."""
  doc = xcif_ext.parse(NULLS_CIF)
  cif_model = iotbx.cif.reader(input_string=NULLS_CIF).model()

  xloop = doc[0].find_loop("_val")
  iblock = list(cif_model.values())[0]

  # String column: raw values should match
  x_notes = xloop.column("_val.note")
  i_notes = list(iblock["_val.note"])
  assert x_notes == i_notes, \
    "notes: xcif=%s vs iotbx=%s" % (x_notes, i_notes)

  # Numeric column: null -> NaN in xcif
  x_vals = xloop.column_as_flex_double("_val.x")
  i_vals_str = iblock["_val.x"]
  assert approx_equal(x_vals[0], 1.5)
  assert math.isnan(x_vals[1]), "Expected NaN for '.'"
  assert math.isnan(x_vals[2]), "Expected NaN for '?'"
  assert approx_equal(x_vals[3], 3.0)

  # Raw string values for nulls
  assert xloop.value(1, 1) == ".", "Expected '.'"
  assert xloop.value(2, 1) == "?", "Expected '?'"
  assert i_vals_str[1] == "."
  assert i_vals_str[2] == "?"


def exercise_quoted_strings():
  """Quoted string values have quotes stripped by both parsers."""
  doc = xcif_ext.parse(QUOTED_CIF)
  cif_model = iotbx.cif.reader(input_string=QUOTED_CIF).model()

  xb = doc[0]
  ib = list(cif_model.values())[0]

  assert xb.find_value("_tag1") == ib["_tag1"], \
    "tag1: xcif='%s' vs iotbx='%s'" % (xb.find_value("_tag1"), ib["_tag1"])
  assert xb.find_value("_tag2") == ib["_tag2"], \
    "tag2: xcif='%s' vs iotbx='%s'" % (xb.find_value("_tag2"), ib["_tag2"])
  assert xb.find_value("_tag3") == ib["_tag3"], \
    "tag3: xcif='%s' vs iotbx='%s'" % (xb.find_value("_tag3"), ib["_tag3"])

  # Verify quotes are stripped
  assert xb.find_value("_tag1") == "single quoted value"
  assert xb.find_value("_tag2") == "double quoted value"
  assert xb.find_value("_tag3") == "unquoted"


def exercise_semicolon_text_field():
  """Semicolon-delimited text fields match between parsers."""
  doc = xcif_ext.parse(SEMICOLON_CIF)
  cif_model = iotbx.cif.reader(input_string=SEMICOLON_CIF).model()

  xb = doc[0]
  ib = list(cif_model.values())[0]

  xtitle = xb.find_value("_struct.title")
  ititle = ib["_struct.title"]

  # Both should contain the multi-line text without the semicolon delimiters
  assert "multi-line" in xtitle, "xcif missing text content: '%s'" % xtitle
  assert "multi-line" in ititle, "iotbx missing text content: '%s'" % ititle

  # Simple tag in the same block
  assert xb.find_value("_other.tag") == "simple"
  assert ib["_other.tag"] == "simple"


def exercise_multi_block():
  """Multi-block CIF files parsed identically."""
  doc = xcif_ext.parse(MULTI_BLOCK_CIF)
  cif_model = iotbx.cif.reader(input_string=MULTI_BLOCK_CIF).model()

  assert len(doc) == len(cif_model), \
    "Block count: xcif=%d vs iotbx=%d" % (len(doc), len(cif_model))

  for i, (bname, iblock) in enumerate(cif_model.items()):
    xblock = doc[i]
    assert xblock.name == bname, \
      "Block %d name: xcif='%s' vs iotbx='%s'" % (i, xblock.name, bname)
    xval = xblock.find_value("_cell.length_a")
    ival = iblock["_cell.length_a"]
    assert xval == ival, \
      "Block %s _cell.length_a: xcif='%s' vs iotbx='%s'" % (bname, xval, ival)


def exercise_find_block_case_insensitive():
  """Case-insensitive block lookup works in both parsers."""
  doc = xcif_ext.parse(MULTI_BLOCK_CIF)
  cif_model = iotbx.cif.reader(input_string=MULTI_BLOCK_CIF).model()

  # xcif: find_block is case-insensitive
  assert doc.find_block("block_a") is not None
  assert doc.find_block("BLOCK_A") is not None
  assert doc.find_block("Block_A") is not None

  # iotbx: dict-like access is case-insensitive
  assert cif_model.get("block_a") is not None
  assert cif_model.get("BLOCK_A") is not None
  assert cif_model.get("Block_A") is not None


def exercise_numeric_conversions():
  """Numeric conversion functions produce correct results."""
  # Standard values
  assert approx_equal(xcif_ext.as_double("10.500"), 10.5)
  assert approx_equal(xcif_ext.as_double("20.000"), 20.0)
  assert approx_equal(xcif_ext.as_double("3.60e-02"), 0.036)

  # SU stripping
  assert approx_equal(xcif_ext.as_double("50.840(10)"), 50.840)
  assert approx_equal(xcif_ext.as_double("42.770(1)"), 42.770)

  # SU extraction
  val, su = xcif_ext.as_double_with_su("50.840(10)")
  assert approx_equal(val, 50.840)
  assert approx_equal(su, 0.010)

  val, su = xcif_ext.as_double_with_su("28.950(2000)")
  assert approx_equal(val, 28.950)
  assert approx_equal(su, 2.000)

  # Null handling
  assert math.isnan(xcif_ext.as_double("."))
  assert math.isnan(xcif_ext.as_double("?"))

  # Integer
  assert xcif_ext.as_int("19") == 19
  assert xcif_ext.as_int("-5") == -5
  assert xcif_ext.as_int("+42") == 42

  # Predicates
  assert xcif_ext.is_null(".")
  assert xcif_ext.is_null("?")
  assert not xcif_ext.is_null("10.0")
  assert xcif_ext.is_unknown(".")
  assert xcif_ext.is_inapplicable("?")


def exercise_flex_column_types():
  """flex.double, flex.int, flex.std_string extraction from loops."""
  doc = xcif_ext.parse(LOOP_CIF)
  loop = doc[0].find_loop("_atom_site")

  # flex.double
  x = loop.column_as_flex_double("_atom_site.Cartn_x")
  assert x.size() == 5
  assert approx_equal(x[0], 27.340)
  assert approx_equal(x[4], 25.112)

  y = loop.column_as_flex_double("_atom_site.Cartn_y")
  assert y.size() == 5
  assert approx_equal(y[0], 24.430)

  # flex.int
  ids = loop.column_as_flex_int("_atom_site.id")
  assert ids.size() == 5
  assert list(ids) == [1, 2, 3, 4, 5]

  # flex.std_string
  syms = loop.column_as_flex_string("_atom_site.type_symbol")
  assert syms.size() == 5
  assert list(syms) == ["N", "C", "C", "O", "C"]

  # Missing column raises KeyError
  try:
    loop.column_as_flex_double("_nonexistent.tag")
    raise AssertionError("Expected KeyError")
  except KeyError:
    pass


def exercise_example_cif():
  """Full comparison on the example.cif fixture file."""
  dist_dir = libtbx.env.dist_path("xcif")
  cif_path = os.path.join(dist_dir, "regression", "example.cif")

  doc = xcif_ext.parse_file(cif_path)
  cif_model = iotbx.cif.reader(file_path=cif_path).model()

  # Block count and name
  assert len(doc) == len(cif_model), \
    "Block count: xcif=%d vs iotbx=%d" % (len(doc), len(cif_model))
  xb = doc[0]
  bname = list(cif_model.keys())[0]
  assert xb.name == bname, \
    "Block name: xcif='%s' vs iotbx='%s'" % (xb.name, bname)
  ib = cif_model[bname]

  # Tag-value pairs
  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma",
              "_symmetry.space_group_name_H-M",
              "_refine.ls_d_res_high", "_refine.ls_R_factor_R_work"]:
    xval = xb.find_value(tag)
    ival = ib[tag]
    assert xval == ival, \
      "%s: xcif='%s' vs iotbx='%s'" % (tag, xval, ival)

  # Null values
  assert xb.find_value("_refine.ls_number_reflns_obs") == "?"
  assert xb.find_value("_refine.B_iso_mean") == "."

  # Atom site loop
  xloop = xb.find_loop("_atom_site")
  iloop = ib.get_loop("_atom_site")
  assert xloop is not None
  assert iloop is not None
  assert xloop.width == iloop.n_columns(), \
    "width: %d vs %d" % (xloop.width, iloop.n_columns())
  assert xloop.length == iloop.n_rows(), \
    "length: %d vs %d" % (xloop.length, iloop.n_rows())

  # Compare all string columns
  for tag in xloop.tags:
    xcol = xloop.column(tag)
    icol = list(iloop[tag])
    assert xcol == icol, \
      "Column %s mismatch: xcif=%s vs iotbx=%s" % (tag, xcol[:3], icol[:3])

  # Compare numeric columns via flex
  x_coords = xloop.column_as_flex_double("_atom_site.Cartn_x")
  i_coords_str = iloop["_atom_site.Cartn_x"]
  for j in range(x_coords.size()):
    assert approx_equal(x_coords[j], float(i_coords_str[j])), \
      "Cartn_x[%d]: %f vs %f" % (j, x_coords[j], float(i_coords_str[j]))


def exercise_save_frames():
  """Save frame access works in xcif."""
  doc = xcif_ext.parse(SAVE_FRAME_CIF)
  block = doc[0]
  assert block.find_value("_dict.title") == "test_dictionary"

  sf = block.find_save_frame("my_category")
  assert sf is not None, "Save frame not found"
  assert sf.find_value("_category.id") == "my_category"
  assert sf.find_value("_category.mandatory_code") == "no"

  assert block.find_save_frame("nonexistent") is None


def exercise_error_on_missing_data_header():
  """Both parsers reject tags outside a data block."""
  # xcif raises RuntimeError
  try:
    xcif_ext.parse("_tag value\n")
    raise AssertionError("xcif: expected exception for missing data_ header")
  except RuntimeError:
    pass

  # iotbx raises CifParserError (subclass of Sorry)
  try:
    iotbx.cif.reader(input_string="_tag value\n")
    raise AssertionError("iotbx: expected exception for missing data_ header")
  except Exception:
    pass


# ---------------------------------------------------------------------------
# Model-level API compatibility: xcif.reader vs iotbx.cif.reader
# Both return model objects accessed through the same dict-like API.
# ---------------------------------------------------------------------------

def exercise_model_simple():
  """xcif.reader and iotbx.cif.reader models match on simple tag-value CIF."""
  m_xcif = xcif.reader(input_string=SIMPLE_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=SIMPLE_CIF).model()

  # Same block names
  assert list(m_xcif.keys()) == list(m_ucif.keys()), \
    "Block names: xcif=%s vs ucif=%s" % (list(m_xcif.keys()), list(m_ucif.keys()))

  bx = m_xcif["compat"]
  bu = m_ucif["compat"]

  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma",
              "_symmetry.space_group_name_H-M"]:
    assert bx[tag] == bu[tag], \
      "Tag %s: xcif='%s' vs ucif='%s'" % (tag, bx[tag], bu[tag])


def exercise_model_loop():
  """xcif.reader and iotbx.cif.reader models match on looped data."""
  m_xcif = xcif.reader(input_string=LOOP_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=LOOP_CIF).model()

  bx = m_xcif["loop_test"]
  bu = m_ucif["loop_test"]

  # Loop dimensions via get_loop
  lx = bx.get_loop("_atom_site")
  lu = bu.get_loop("_atom_site")
  assert lx is not None
  assert lu is not None
  assert lx.n_rows() == lu.n_rows(), \
    "n_rows: xcif=%d vs ucif=%d" % (lx.n_rows(), lu.n_rows())
  assert lx.n_columns() == lu.n_columns(), \
    "n_columns: xcif=%d vs ucif=%d" % (lx.n_columns(), lu.n_columns())

  # Column access via block[tag] -> flex.std_string
  for tag in ["_atom_site.id", "_atom_site.type_symbol",
              "_atom_site.label_comp_id",
              "_atom_site.Cartn_x", "_atom_site.Cartn_y",
              "_atom_site.Cartn_z", "_atom_site.occupancy",
              "_atom_site.B_iso_or_equiv"]:
    col_x = list(bx[tag])
    col_u = list(bu[tag])
    assert col_x == col_u, \
      "Column %s: xcif=%s vs ucif=%s" % (tag, col_x[:3], col_u[:3])


def exercise_model_multi_block():
  """xcif.reader and iotbx.cif.reader models match on multi-block CIF."""
  m_xcif = xcif.reader(input_string=MULTI_BLOCK_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=MULTI_BLOCK_CIF).model()

  assert list(m_xcif.keys()) == list(m_ucif.keys())
  assert len(m_xcif) == len(m_ucif)

  for name in m_ucif.keys():
    bx = m_xcif[name]
    bu = m_ucif[name]
    assert bx["_cell.length_a"] == bu["_cell.length_a"], \
      "Block %s _cell.length_a: xcif='%s' vs ucif='%s'" % (
        name, bx["_cell.length_a"], bu["_cell.length_a"])
    assert bx["_entry.id"] == bu["_entry.id"]


def exercise_model_nulls():
  """xcif.reader and iotbx.cif.reader models match on null values."""
  m_xcif = xcif.reader(input_string=NULLS_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=NULLS_CIF).model()

  bx = m_xcif["nulls"]
  bu = m_ucif["nulls"]

  for tag in ["_val.id", "_val.x", "_val.note"]:
    col_x = list(bx[tag])
    col_u = list(bu[tag])
    assert col_x == col_u, \
      "Column %s: xcif=%s vs ucif=%s" % (tag, col_x, col_u)


def exercise_model_quoted():
  """xcif.reader and iotbx.cif.reader models match on quoted strings."""
  m_xcif = xcif.reader(input_string=QUOTED_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=QUOTED_CIF).model()

  bx = m_xcif["quoted"]
  bu = m_ucif["quoted"]

  for tag in ["_tag1", "_tag2", "_tag3"]:
    assert bx[tag] == bu[tag], \
      "Tag %s: xcif='%s' vs ucif='%s'" % (tag, bx[tag], bu[tag])


def exercise_model_case_insensitive():
  """Both models support case-insensitive block and tag lookup."""
  m_xcif = xcif.reader(input_string=SIMPLE_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=SIMPLE_CIF).model()

  # Case-insensitive block lookup
  assert "COMPAT" in m_xcif
  assert "COMPAT" in m_ucif
  assert "compat" in m_xcif
  assert "compat" in m_ucif

  # Case-insensitive tag lookup
  bx = m_xcif["COMPAT"]
  bu = m_ucif["COMPAT"]
  assert bx["_CELL.LENGTH_A"] == bu["_CELL.LENGTH_A"]
  assert bx["_cell.length_a"] == bu["_cell.length_a"]


def exercise_model_example_cif():
  """Full model comparison on example.cif fixture."""
  dist_dir = libtbx.env.dist_path("xcif")
  cif_path = os.path.join(dist_dir, "regression", "example.cif")

  m_xcif = xcif.reader(file_path=cif_path).model()
  m_ucif = iotbx.cif.reader(file_path=cif_path).model()

  # Same block names
  assert list(m_xcif.keys()) == list(m_ucif.keys())

  bx = m_xcif["1UBQ"]
  bu = m_ucif["1UBQ"]

  # Tag-value pairs match
  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma",
              "_symmetry.space_group_name_H-M",
              "_refine.ls_d_res_high", "_refine.ls_R_factor_R_work",
              "_refine.ls_R_factor_R_free",
              "_refine.ls_number_reflns_obs", "_refine.B_iso_mean"]:
    assert bx[tag] == bu[tag], \
      "Tag %s: xcif='%s' vs ucif='%s'" % (tag, bx[tag], bu[tag])

  # Loop dimensions match
  lx = bx.get_loop("_atom_site")
  lu = bu.get_loop("_atom_site")
  assert lx.n_rows() == lu.n_rows()
  assert lx.n_columns() == lu.n_columns()

  # All loop columns match cell-by-cell
  for tag in lx.keys():
    col_x = list(bx[tag])
    col_u = list(bu[tag])
    assert col_x == col_u, \
      "Column %s: xcif=%s vs ucif=%s" % (tag, col_x[:3], col_u[:3])

  # Reflns loop
  lx_r = bx.get_loop("_reflns")
  lu_r = bu.get_loop("_reflns")
  assert lx_r is not None
  assert lu_r is not None
  assert lx_r.n_rows() == lu_r.n_rows()
  for tag in lx_r.keys():
    col_x = list(bx[tag])
    col_u = list(bu[tag])
    assert col_x == col_u, \
      "Reflns column %s: xcif=%s vs ucif=%s" % (tag, col_x, col_u)


# ---------------------------------------------------------------------------
# Identical show() output: xcif.reader vs iotbx.cif.reader
# ---------------------------------------------------------------------------

try:
  from six.moves import StringIO
except ImportError:
  from io import StringIO


def _show_str(model):
  sio = StringIO()
  model.show(out=sio)
  return sio.getvalue()


def exercise_output_simple():
  """show() output matches for simple tag-value CIF."""
  m_xcif = xcif.reader(input_string=SIMPLE_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=SIMPLE_CIF).model()
  sx = _show_str(m_xcif)
  su = _show_str(m_ucif)
  assert sx == su, "Simple CIF output differs:\nxcif:\n%s\nucif:\n%s" % (sx, su)


def exercise_output_loop():
  """show() output matches for looped data."""
  m_xcif = xcif.reader(input_string=LOOP_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=LOOP_CIF).model()
  sx = _show_str(m_xcif)
  su = _show_str(m_ucif)
  assert sx == su, "Loop CIF output differs:\nxcif:\n%s\nucif:\n%s" % (sx, su)


def exercise_output_multi_block():
  """show() output matches for multi-block CIF."""
  m_xcif = xcif.reader(input_string=MULTI_BLOCK_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=MULTI_BLOCK_CIF).model()
  sx = _show_str(m_xcif)
  su = _show_str(m_ucif)
  assert sx == su, \
    "Multi-block output differs:\nxcif:\n%s\nucif:\n%s" % (sx, su)


def exercise_output_quoted():
  """show() output matches for quoted strings."""
  m_xcif = xcif.reader(input_string=QUOTED_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=QUOTED_CIF).model()
  sx = _show_str(m_xcif)
  su = _show_str(m_ucif)
  assert sx == su, \
    "Quoted string output differs:\nxcif:\n%s\nucif:\n%s" % (sx, su)


def exercise_output_nulls():
  """show() output matches for null values."""
  m_xcif = xcif.reader(input_string=NULLS_CIF).model()
  m_ucif = iotbx.cif.reader(input_string=NULLS_CIF).model()
  sx = _show_str(m_xcif)
  su = _show_str(m_ucif)
  assert sx == su, \
    "Nulls CIF output differs:\nxcif:\n%s\nucif:\n%s" % (sx, su)


def exercise_output_example_cif():
  """show() output matches for the full example.cif fixture."""
  dist_dir = libtbx.env.dist_path("xcif")
  cif_path = os.path.join(dist_dir, "regression", "example.cif")
  m_xcif = xcif.reader(file_path=cif_path).model()
  m_ucif = iotbx.cif.reader(file_path=cif_path).model()
  sx = _show_str(m_xcif)
  su = _show_str(m_ucif)
  assert sx == su, \
    "example.cif output differs:\nxcif:\n%s\nucif:\n%s" % (
      sx[:500], su[:500])


def exercise_loop_category():
  """Edge cases for the _loop_category() helper used by _loop_adapter."""
  from xcif import _loop_category

  # Empty list -> empty string
  assert _loop_category([]) == "", \
    "Expected '' for empty list, got '%s'" % _loop_category([])

  # Single tag returned as-is (no trimming applied)
  assert _loop_category(["_atom_site.id"]) == "_atom_site.id", \
    "Single tag with dot: %s" % _loop_category(["_atom_site.id"])
  assert _loop_category(["_single"]) == "_single", \
    "Single tag no dot: %s" % _loop_category(["_single"])

  # Multiple tags sharing a dot-separated category
  r = _loop_category(["_atom_site.id", "_atom_site.x", "_atom_site.z"])
  assert r == "_atom_site", \
    "Dot-prefix: expected '_atom_site', got '%s'" % r

  # Tags that share only the category name (diverge right after the dot)
  # "_cell.length_a" vs "_cell.angle_alpha" → common prefix = "_cell." → "_cell"
  r = _loop_category(["_cell.length_a", "_cell.angle_alpha"])
  assert r == "_cell", "Expected '_cell', got '%s'" % r

  # Multiple tags sharing only an underscore delimiter (no dot)
  # "_refln_h", "_refln_k" → common prefix "_refln_" (ends with '_', kept)
  r = _loop_category(["_refln_h", "_refln_k", "_refln_l"])
  assert r.startswith("_refln"), \
    "Underscore-prefix: expected '_refln...', got '%s'" % r


if __name__ == "__main__":
  # C++ binding level tests
  exercise_block_names()
  exercise_tag_values()
  exercise_tag_case_insensitive()
  exercise_loop_columns()
  exercise_loop_dimensions()
  exercise_null_handling()
  exercise_quoted_strings()
  exercise_semicolon_text_field()
  exercise_multi_block()
  exercise_find_block_case_insensitive()
  exercise_numeric_conversions()
  exercise_flex_column_types()
  exercise_example_cif()
  exercise_save_frames()
  exercise_error_on_missing_data_header()
  # Model-level: xcif.reader vs iotbx.cif.reader
  exercise_model_simple()
  exercise_model_loop()
  exercise_model_multi_block()
  exercise_model_nulls()
  exercise_model_quoted()
  exercise_model_case_insensitive()
  exercise_model_example_cif()
  # Identical output: show() comparison
  exercise_output_simple()
  exercise_output_loop()
  exercise_output_multi_block()
  exercise_output_quoted()
  exercise_output_nulls()
  exercise_output_example_cif()
  # Internal helper
  exercise_loop_category()
  print(format_cpu_times())
  print("OK")
