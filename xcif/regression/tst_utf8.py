from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_utf8.py
#
# Tests that xcif handles non-ASCII (UTF-8) content gracefully.
# CIF 1.1 restricts characters to ASCII printable (32-126), but
# real-world files often contain UTF-8 in comments, quoted strings,
# and semicolon text fields.  xcif should parse these without errors.

from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex  # registers shared<double> converter


def exercise_utf8_in_comments():
  """UTF-8 characters in comments are skipped without error."""
  import xcif_ext

  cif_text = (
    u"# This comment has an em-dash \u2014 and a degree sign \u00b0\n"
    u"data_utf8_comment\n"
    u"_cell.length_a 10.0\n"
    u"# Greek: \u03b1 \u03b2 \u03b3\n"
    u"_cell.length_b 20.0\n"
  ).encode("utf-8").decode("utf-8")

  doc = xcif_ext.parse(cif_text)
  assert len(doc) == 1
  block = doc[0]
  assert block.name == "utf8_comment"
  assert block.find_value("_cell.length_a") == "10.0"
  assert block.find_value("_cell.length_b") == "20.0"


def exercise_utf8_in_quoted_string():
  """UTF-8 characters in quoted strings are preserved."""
  import xcif_ext

  cif_text = (
    u"data_utf8_quoted\n"
    u"_author.name 'M\u00fcller'\n"
    u"_journal.name \"Zeitschrift f\u00fcr Kristallographie\"\n"
  ).encode("utf-8").decode("utf-8")

  doc = xcif_ext.parse(cif_text)
  block = doc[0]
  assert block.find_value("_author.name") == u"M\u00fcller"
  assert u"f\u00fcr" in block.find_value("_journal.name")


def exercise_utf8_in_semicolon_field():
  """UTF-8 characters in semicolon text fields are preserved."""
  import xcif_ext

  cif_text = (
    u"data_utf8_semi\n"
    u"_struct.title\n"
    u";Crystal structure of the \u03b1-polymorph.\n"
    u"Measured at 20 \u00b0C.\n"
    u";\n"
    u"_other.tag ok\n"
  ).encode("utf-8").decode("utf-8")

  doc = xcif_ext.parse(cif_text)
  block = doc[0]

  title = block.find_value("_struct.title")
  assert u"\u03b1" in title, "Missing alpha: '%s'" % title
  assert u"\u00b0" in title, "Missing degree sign: '%s'" % title
  assert block.find_value("_other.tag") == "ok"


def exercise_utf8_in_loop_values():
  """UTF-8 characters in loop cell values are preserved."""
  import xcif_ext

  cif_text = (
    u"data_utf8_loop\n"
    u"loop_\n"
    u"_author.id\n"
    u"_author.name\n"
    u"1 'Sch\u00f6n'\n"
    u"2 'Bj\u00f6rk'\n"
    u"3 'Garc\u00eda'\n"
  ).encode("utf-8").decode("utf-8")

  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_author")
  assert loop is not None
  assert loop.length == 3

  names = loop.column("_author.name")
  assert names[0] == u"Sch\u00f6n", "Got '%s'" % names[0]
  assert names[1] == u"Bj\u00f6rk", "Got '%s'" % names[1]
  assert names[2] == u"Garc\u00eda", "Got '%s'" % names[2]


def exercise_utf8_block_data_integrity():
  """Numeric data after UTF-8 comments parses correctly."""
  import xcif_ext

  cif_text = (
    u"# Crystallographic data for compound \u03b1-Fe\u2082O\u2083\n"
    u"data_alpha_fe2o3\n"
    u"# Unit cell (\u00c5ngstr\u00f6ms)\n"
    u"_cell.length_a 5.038\n"
    u"_cell.length_b 5.038\n"
    u"_cell.length_c 13.772\n"
    u"loop_\n"
    u"_atom.id\n"
    u"_atom.x\n"
    u"_atom.y\n"
    u"1 0.000 0.000\n"
    u"2 0.333 0.333\n"
    u"3 0.667 0.667\n"
  ).encode("utf-8").decode("utf-8")

  doc = xcif_ext.parse(cif_text)
  block = doc[0]
  assert block.name == "alpha_fe2o3"
  assert approx_equal(xcif_ext.as_double(block.find_value("_cell.length_a")),
                       5.038)
  assert approx_equal(xcif_ext.as_double(block.find_value("_cell.length_c")),
                       13.772)

  loop = block.find_loop("_atom")
  x = loop.column_as_flex_double("_atom.x")
  assert approx_equal(list(x), [0.0, 0.333, 0.667])


def exercise_pure_ascii_unaffected():
  """Pure ASCII CIF still works correctly (sanity check)."""
  import xcif_ext

  cif_text = (
    "data_ascii\n"
    "_tag.a 1.0\n"
    "_tag.b 'hello world'\n"
    "loop_\n"
    "_val.x\n"
    "1.0\n"
    "2.0\n"
    "3.0\n"
  )
  doc = xcif_ext.parse(cif_text)
  assert doc[0].find_value("_tag.a") == "1.0"
  assert doc[0].find_value("_tag.b") == "hello world"
  loop = doc[0].find_loop("_val")
  assert loop.length == 3


if __name__ == "__main__":
  exercise_utf8_in_comments()
  exercise_utf8_in_quoted_string()
  exercise_utf8_in_semicolon_field()
  exercise_utf8_in_loop_values()
  exercise_utf8_block_data_integrity()
  exercise_pure_ascii_unaffected()
  print(format_cpu_times())
  print("OK")
