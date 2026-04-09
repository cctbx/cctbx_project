from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_bindings.py
#
# Core Python binding tests for the xcif_ext Boost.Python module.
# Tests the complete Python API surface: Document, Block, Loop,
# numeric functions, flex column extraction, and error handling.
#
# These tests define the DESIRED API for Step 10 implementation.
# They will fail with ImportError until xcif_ext is built.

import os
import math
import libtbx.load_env
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal

def exercise_module_attributes():
  """Verify xcif_ext exposes the expected top-level symbols."""
  import xcif_ext
  for name in ('parse', 'parse_file',
               'as_double', 'as_int', 'as_double_with_su',
               'is_null', 'is_unknown', 'is_inapplicable'):
    assert hasattr(xcif_ext, name), "Missing attribute: %s" % name

def exercise_parse_string():
  """Parse a CIF string, navigate Document -> Block -> Loop."""
  import xcif_ext
  cif_text = (
    "data_test\n"
    "_cell.length_a 10.5\n"
    "_cell.length_b 20.0\n"
    "loop_\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.Cartn_x\n"
    "1 C 1.5\n"
    "2 N 3.5\n"
    "3 O 5.5\n"
  )
  doc = xcif_ext.parse(cif_text)
  assert len(doc) == 1, "Expected 1 block, got %d" % len(doc)

  block = doc[0]
  assert block.name == "test", "Expected 'test', got '%s'" % block.name

  # find_block: case-insensitive
  assert doc.find_block("test") is not None
  assert doc.find_block("TEST") is not None
  assert doc.find_block("nonexistent") is None

def exercise_tag_values():
  """Tag-value pair access and case-insensitive lookup."""
  import xcif_ext
  cif_text = (
    "data_tv\n"
    "_cell.length_a 10.5\n"
    "_symmetry.space_group_name_H-M 'P 21 21 21'\n"
  )
  doc = xcif_ext.parse(cif_text)
  block = doc[0]

  assert block.has_tag("_cell.length_a")
  assert block.has_tag("_CELL.LENGTH_A")
  assert not block.has_tag("_nonexistent")

  val = block.find_value("_cell.length_a")
  assert val == "10.5", "Expected '10.5', got '%s'" % val

  sg = block.find_value("_symmetry.space_group_name_H-M")
  assert sg == "P 21 21 21", "Expected 'P 21 21 21', got '%s'" % sg

  # Missing tag returns None
  assert block.find_value("_nonexistent") is None

def exercise_loop_access():
  """Loop navigation, dimensions, cell access, and column extraction."""
  import xcif_ext
  cif_text = (
    "data_loop\n"
    "loop_\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "1 C 1.5 2.5\n"
    "2 N 3.5 4.5\n"
    "3 O 5.5 6.5\n"
  )
  doc = xcif_ext.parse(cif_text)
  block = doc[0]

  # find_loop by category prefix
  loop = block.find_loop("_atom_site")
  assert loop is not None

  # find_loop by full tag name
  assert block.find_loop("_atom_site.id") is not None
  assert block.find_loop("_nonexistent") is None

  # Dimensions
  assert loop.width == 4, "Expected width 4, got %d" % loop.width
  assert loop.length == 3, "Expected length 3, got %d" % loop.length

  # Tags list
  tags = loop.tags
  assert len(tags) == 4
  assert tags[0] == "_atom_site.id"

  # has_tag / column_index
  assert loop.has_tag("_atom_site.Cartn_x")
  assert not loop.has_tag("_nonexistent")
  ci = loop.column_index("_atom_site.Cartn_x")
  assert ci == 2, "Expected column_index 2, got %d" % ci

  # Missing column_index returns -1
  assert loop.column_index("_nonexistent") == -1

  # Cell value access
  assert loop.value(0, 2) == "1.5"
  assert loop.value(2, 1) == "O"

  # Column extraction as Python list of strings
  col = loop.column("_atom_site.type_symbol")
  assert col == ["C", "N", "O"], "Got %s" % col

  # loops property on block
  loops = block.loops
  assert len(loops) == 1

def exercise_flex_columns():
  """Direct flex array extraction from loops (no intermediate Python list)."""
  import xcif_ext
  from scitbx.array_family import flex

  cif_text = (
    "data_flex\n"
    "loop_\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.Cartn_x\n"
    "1 C 1.5\n"
    "2 N 3.5\n"
    "3 O 5.5\n"
  )
  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_atom_site")

  # flex.double
  x = loop.column_as_flex_double("_atom_site.Cartn_x")
  assert x.size() == 3
  assert approx_equal(x[0], 1.5)
  assert approx_equal(x[1], 3.5)
  assert approx_equal(x[2], 5.5)

  # flex.int
  ids = loop.column_as_flex_int("_atom_site.id")
  assert ids.size() == 3
  assert ids[0] == 1 and ids[1] == 2 and ids[2] == 3

  # flex.std_string
  syms = loop.column_as_flex_string("_atom_site.type_symbol")
  assert syms.size() == 3
  assert syms[0] == "C" and syms[1] == "N" and syms[2] == "O"

  # flex.double with null values -> NaN
  cif_null = (
    "data_nulls\n"
    "loop_\n"
    "_val.x\n"
    "1.0\n"
    ".\n"
    "?\n"
    "3.0\n"
  )
  doc2 = xcif_ext.parse(cif_null)
  loop2 = doc2[0].find_loop("_val")
  x2 = loop2.column_as_flex_double("_val.x")
  assert x2.size() == 4
  assert approx_equal(x2[0], 1.0)
  assert math.isnan(x2[1]), "Expected NaN for '.'"
  assert math.isnan(x2[2]), "Expected NaN for '?'"
  assert approx_equal(x2[3], 3.0)

  # column_as_flex_int raises ValueError when a cell contains '.' or '?'
  cif_nullints = (
    "data_nullints\n"
    "loop_\n"
    "_val.n\n"
    "1\n"
    ".\n"
    "3\n"
  )
  doc3 = xcif_ext.parse(cif_nullints)
  loop3 = doc3[0].find_loop("_val")
  try:
    loop3.column_as_flex_int("_val.n")
    raise AssertionError("Expected ValueError for '.' in int column")
  except ValueError:
    pass

def exercise_numeric():
  """Numeric conversion free functions."""
  import xcif_ext

  # as_double
  assert approx_equal(xcif_ext.as_double("50.840"), 50.840)
  assert approx_equal(xcif_ext.as_double("3.60e-02"), 0.036)
  assert approx_equal(xcif_ext.as_double("50.840(10)"), 50.840)
  assert math.isnan(xcif_ext.as_double("."))
  assert math.isnan(xcif_ext.as_double("?"))

  # as_int
  assert xcif_ext.as_int("19") == 19
  assert xcif_ext.as_int("-5") == -5
  assert xcif_ext.as_int("+42") == 42
  try:
    xcif_ext.as_int("abc")
    raise AssertionError("Expected ValueError for non-integer")
  except (ValueError, RuntimeError):
    pass

  # as_double_with_su
  val, su = xcif_ext.as_double_with_su("50.840(10)")
  assert approx_equal(val, 50.840)
  assert approx_equal(su, 0.010)
  val, su = xcif_ext.as_double_with_su("42.770(1)")
  assert approx_equal(val, 42.770)
  assert approx_equal(su, 0.001)
  val, su = xcif_ext.as_double_with_su("10.0")
  assert approx_equal(val, 10.0)
  assert approx_equal(su, 0.0)

  # Predicates
  assert xcif_ext.is_null(".")
  assert xcif_ext.is_null("?")
  assert not xcif_ext.is_null("10.0")
  assert xcif_ext.is_unknown(".")
  assert not xcif_ext.is_unknown("?")
  assert xcif_ext.is_inapplicable("?")
  assert not xcif_ext.is_inapplicable(".")

def exercise_parse_file():
  """File-based parsing (mmap path) with example.cif."""
  import xcif_ext
  dist_dir = libtbx.env.dist_path("xcif")
  cif_file = os.path.join(dist_dir, "regression", "example.cif")

  doc = xcif_ext.parse_file(cif_file)
  assert len(doc) == 1
  block = doc[0]
  assert block.name.upper() == "1UBQ"

  val = block.find_value("_cell.length_a")
  assert val == "50.840(10)", "Expected '50.840(10)', got '%s'" % val

  loop = block.find_loop("_atom_site")
  assert loop is not None
  assert loop.width == 13
  assert loop.length == 10

def exercise_error_handling():
  """Parse errors raise Python exceptions with location info."""
  import xcif_ext

  # Tag outside a data block
  try:
    xcif_ext.parse("_tag value")
    raise AssertionError("Expected exception for missing data_ header")
  except RuntimeError as e:
    assert "line" in str(e).lower() or "1" in str(e)

  # Loop with mismatched value count
  try:
    xcif_ext.parse("data_x\nloop_\n_a\n_b\nv1\n")
    raise AssertionError("Expected exception for mismatched loop values")
  except RuntimeError:
    pass

  # Index out of range on Document
  doc = xcif_ext.parse("data_a\n")
  try:
    doc[99]
    raise AssertionError("Expected IndexError")
  except IndexError:
    pass

def exercise_save_frames():
  """Save frame access within a block."""
  import xcif_ext
  cif_text = (
    "data_dict\n"
    "save_my_frame\n"
    "_tag1 value1\n"
    "save_\n"
  )
  doc = xcif_ext.parse(cif_text)
  block = doc[0]
  sf = block.find_save_frame("my_frame")
  assert sf is not None
  assert sf.find_value("_tag1") == "value1"
  assert block.find_save_frame("nonexistent") is None

def exercise_multi_block():
  """Multiple data blocks in one document."""
  import xcif_ext
  cif_text = "data_a\n_x 1\ndata_b\n_y 2\ndata_c\n_z 3\n"
  doc = xcif_ext.parse(cif_text)
  assert len(doc) == 3
  assert doc[0].name == "a"
  assert doc[1].name == "b"
  assert doc[2].name == "c"
  assert doc.find_block("b").find_value("_y") == "2"

if __name__ == "__main__":
  exercise_module_attributes()
  exercise_parse_string()
  exercise_tag_values()
  exercise_loop_access()
  exercise_flex_columns()
  exercise_numeric()
  exercise_parse_file()
  exercise_error_handling()
  exercise_save_frames()
  exercise_multi_block()
  print(format_cpu_times())
  print("OK")
