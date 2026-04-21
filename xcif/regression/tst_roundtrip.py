from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_roundtrip.py
#
# Step 12: True round-trip tests for xcif.reader and adapter show().
# Parse CIF -> xcif.reader model -> show() -> re-parse -> compare.

import os
import xcif
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import libtbx.load_env

try:
  from six.moves import StringIO
except ImportError:
  from io import StringIO


def _model_to_string(model):
  sio = StringIO()
  model.show(out=sio)
  return sio.getvalue()


def exercise_reader_input_string():
  """xcif.reader(input_string=...) returns a usable model."""
  cif_text = (
    "data_test\n"
    "_cell.length_a 10.5\n"
    "_cell.length_b 20.0\n"
  )
  r = xcif.reader(input_string=cif_text)
  m = r.model()
  assert len(m) == 1
  assert "test" in m
  block = m["test"]
  assert block["_cell.length_a"] == "10.5"
  assert block["_cell.length_b"] == "20.0"


def exercise_reader_file_object():
  """xcif.reader(file_object=...) works."""
  cif_text = (
    "data_fo\n"
    "_tag.x 42\n"
  )
  r = xcif.reader(file_object=StringIO(cif_text))
  m = r.model()
  assert "fo" in m
  assert m["fo"]["_tag.x"] == "42"


def exercise_reader_file_path():
  """xcif.reader(file_path=...) works on example.cif."""
  dist_dir = libtbx.env.dist_path("xcif")
  cif_path = os.path.join(dist_dir, "regression", "example.cif")
  r = xcif.reader(file_path=cif_path)
  m = r.model()
  assert len(m) == 1
  block = m["1UBQ"]
  assert approx_equal(float(block["_cell.angle_alpha"]), 90.0)


def exercise_roundtrip_tag_values():
  """Tag-value pairs survive show() -> re-parse round-trip."""
  cif_text = (
    "data_tv\n"
    "_cell.length_a 25.000\n"
    "_cell.length_b 30.000\n"
    "_cell.length_c 35.000\n"
    "_symmetry.space_group_name_H-M 'P 21 21 21'\n"
  )
  m1 = xcif.reader(input_string=cif_text).model()
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()

  b1 = m1["tv"]
  b2 = m2["tv"]
  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_symmetry.space_group_name_H-M"]:
    assert b1[tag] == b2[tag], \
      "Tag %s: '%s' != '%s'" % (tag, b1[tag], b2[tag])


def exercise_roundtrip_loop():
  """Loop data survives show() -> re-parse round-trip."""
  cif_text = (
    "data_lp\n"
    "loop_\n"
    "_atom.id\n"
    "_atom.symbol\n"
    "_atom.x\n"
    "_atom.y\n"
    "1 C 1.000 2.000\n"
    "2 N 3.000 4.000\n"
    "3 O 5.000 6.000\n"
  )
  m1 = xcif.reader(input_string=cif_text).model()
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()

  b1 = m1["lp"]
  b2 = m2["lp"]

  lp1 = b1.get_loop("_atom")
  lp2 = b2.get_loop("_atom")
  assert lp1 is not None
  assert lp2 is not None
  assert lp1.n_rows() == lp2.n_rows()
  assert lp1.n_columns() == lp2.n_columns()

  for tag in ["_atom.id", "_atom.symbol", "_atom.x", "_atom.y"]:
    col1 = list(b1[tag])
    col2 = list(b2[tag])
    assert col1 == col2, "Column %s: %s != %s" % (tag, col1, col2)


def exercise_roundtrip_quoted_strings():
  """Quoted strings survive round-trip (quotes stripped, re-added)."""
  cif_text = (
    "data_qs\n"
    "_tag1 'value with spaces'\n"
    "_tag2 \"double quoted\"\n"
    "_tag3 simple\n"
    "_tag4 'P 1 21 1'\n"
  )
  m1 = xcif.reader(input_string=cif_text).model()
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()

  b1 = m1["qs"]
  b2 = m2["qs"]
  assert b1["_tag1"] == b2["_tag1"] == "value with spaces"
  assert b1["_tag2"] == b2["_tag2"] == "double quoted"
  assert b1["_tag3"] == b2["_tag3"] == "simple"
  assert b1["_tag4"] == b2["_tag4"] == "P 1 21 1"


def exercise_roundtrip_semicolon_field():
  """Semicolon text fields survive round-trip.

  Like iotbx, show() adds a leading newline to semicolon field content.
  First pass: value changes (leading \\n added). Second pass: stable.
  This matches iotbx.cif.model round-trip behavior exactly.
  """
  cif_text = (
    "data_sf\n"
    "_struct.title\n"
    ";First line of text.\n"
    "Second line of text.\n"
    ";\n"
    "_other.tag value\n"
  )
  m1 = xcif.reader(input_string=cif_text).model()
  b1 = m1["sf"]
  assert "First line" in b1["_struct.title"]
  assert "Second line" in b1["_struct.title"]

  # First round-trip: show() adds leading \n (same as iotbx behavior)
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()
  b2 = m2["sf"]

  # Second round-trip: should now be stable
  text_out2 = _model_to_string(m2)
  m3 = xcif.reader(input_string=text_out2).model()
  b3 = m3["sf"]
  assert b2["_struct.title"] == b3["_struct.title"], \
    "Semicolon field not stable after second round-trip:\n'%s'\nvs\n'%s'" % (
      b2["_struct.title"], b3["_struct.title"])
  assert b2["_other.tag"] == b3["_other.tag"]
  assert b1["_other.tag"] == b2["_other.tag"]


def exercise_roundtrip_null_values():
  """Null/unknown/inapplicable values survive round-trip."""
  cif_text = (
    "data_nv\n"
    "_tag1 .\n"
    "_tag2 ?\n"
    "loop_\n"
    "_val.id\n"
    "_val.x\n"
    "1 .\n"
    "2 ?\n"
    "3 1.5\n"
  )
  m1 = xcif.reader(input_string=cif_text).model()
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()

  b1 = m1["nv"]
  b2 = m2["nv"]
  assert b1["_tag1"] == b2["_tag1"] == "."
  assert b1["_tag2"] == b2["_tag2"] == "?"

  col1 = list(b1["_val.x"])
  col2 = list(b2["_val.x"])
  assert col1 == col2 == [".", "?", "1.5"]


def exercise_roundtrip_multi_block():
  """Multi-block document survives round-trip, block order preserved."""
  cif_text = (
    "data_alpha\n"
    "_tag.a 1\n"
    "data_beta\n"
    "_tag.b 2\n"
    "data_gamma\n"
    "_tag.c 3\n"
  )
  m1 = xcif.reader(input_string=cif_text).model()
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()

  names1 = list(m1.keys())
  names2 = list(m2.keys())
  assert names1 == names2, "%s != %s" % (names1, names2)
  assert m2["alpha"]["_tag.a"] == "1"
  assert m2["beta"]["_tag.b"] == "2"
  assert m2["gamma"]["_tag.c"] == "3"


def exercise_roundtrip_case_insensitive():
  """Case-insensitive block/tag lookup works through the adapter."""
  cif_text = (
    "data_MixedCase\n"
    "_Cell.Length_A 10.0\n"
  )
  m = xcif.reader(input_string=cif_text).model()
  assert "mixedcase" in m
  assert "MIXEDCASE" in m
  assert "MixedCase" in m
  block = m["mixedcase"]
  assert block["_cell.length_a"] == "10.0"
  assert block["_CELL.LENGTH_A"] == "10.0"


def exercise_roundtrip_example_cif():
  """Full example.cif round-trip: parse -> show -> re-parse -> compare."""
  dist_dir = libtbx.env.dist_path("xcif")
  cif_path = os.path.join(dist_dir, "regression", "example.cif")

  m1 = xcif.reader(file_path=cif_path).model()
  text_out = _model_to_string(m1)
  m2 = xcif.reader(input_string=text_out).model()

  # Block names match
  assert list(m1.keys()) == list(m2.keys())

  b1 = m1["1UBQ"]
  b2 = m2["1UBQ"]

  # Tag-value pairs match
  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma",
              "_symmetry.space_group_name_H-M"]:
    assert b1[tag] == b2[tag], \
      "Tag %s: '%s' != '%s'" % (tag, b1[tag], b2[tag])

  # Semicolon text field: show() adds leading \n on first pass,
  # so compare content rather than exact string
  assert "ubiquitin" in b2["_struct.title"].lower()
  # Second round-trip should be stable
  text_out2 = _model_to_string(m2)
  m3 = xcif.reader(input_string=text_out2).model()
  assert m2["1UBQ"]["_struct.title"] == m3["1UBQ"]["_struct.title"]

  # Atom site loop match
  lp1 = b1.get_loop("_atom_site")
  lp2 = b2.get_loop("_atom_site")
  assert lp1 is not None
  assert lp2 is not None
  assert lp1.n_rows() == lp2.n_rows()
  assert lp1.n_columns() == lp2.n_columns()

  for tag in lp1.keys():
    col1 = list(b1[tag])
    col2 = list(b2[tag])
    assert col1 == col2, "Column %s mismatch" % tag


def exercise_block_loops_property():
  """block.loops property returns a dict of category -> loop."""
  cif_text = (
    "data_lps\n"
    "_tag.x 1\n"
    "loop_\n"
    "_atom.id\n"
    "_atom.symbol\n"
    "1 C\n"
    "2 N\n"
    "loop_\n"
    "_reflns.id\n"
    "_reflns.d_max\n"
    "1 50.0\n"
  )
  m = xcif.reader(input_string=cif_text).model()
  b = m["lps"]
  loops = b.loops
  assert len(loops) >= 2


def exercise_loop_adapter_methods():
  """Loop adapter size/n_rows/n_columns/name work."""
  cif_text = (
    "data_la\n"
    "loop_\n"
    "_atom.id\n"
    "_atom.symbol\n"
    "_atom.x\n"
    "1 C 1.0\n"
    "2 N 2.0\n"
  )
  m = xcif.reader(input_string=cif_text).model()
  lp = m["la"].get_loop("_atom")
  assert lp is not None
  assert lp.n_rows() == 2
  assert lp.n_columns() == 3
  assert lp.size() == 2
  name = lp.name()
  assert name is not None


def exercise_model_show_str():
  """str() on model and block produces valid CIF text."""
  cif_text = (
    "data_s\n"
    "_tag.a 1\n"
  )
  m = xcif.reader(input_string=cif_text).model()
  s = str(m)
  assert "data_s" in s
  assert "_tag.a" in s

  bs = str(m["s"])
  assert "_tag.a" in bs


def exercise_error_count():
  """reader.error_count() and show_errors() work."""
  r = xcif.reader(input_string="data_e\n_t 1\n")
  assert r.error_count() == 0
  r.show_errors()  # should not raise


if __name__ == "__main__":
  exercise_reader_input_string()
  exercise_reader_file_object()
  exercise_reader_file_path()
  exercise_roundtrip_tag_values()
  exercise_roundtrip_loop()
  exercise_roundtrip_quoted_strings()
  exercise_roundtrip_semicolon_field()
  exercise_roundtrip_null_values()
  exercise_roundtrip_multi_block()
  exercise_roundtrip_case_insensitive()
  exercise_roundtrip_example_cif()
  exercise_block_loops_property()
  exercise_loop_adapter_methods()
  exercise_model_show_str()
  exercise_error_count()
  print(format_cpu_times())
  print("OK")
