"""Tests for xcif_ext.parse(..., strict=False) — synthesized global_ block
for content before the first data_ block.

Motivating case: the cctbx monomer library files start with pair items
(e.g. _lib_update 15/04/05) before the first data_ block. ucif silently
accepts this in non-strict mode; xcif must also, attaching the leading
content to an implicit global_ block.
"""
from __future__ import absolute_import, division, print_function

import xcif_ext


def _expect_parse_error(input_string, **kwargs):
  try:
    xcif_ext.parse(input_string, **kwargs)
  except (RuntimeError, ValueError):
    return
  raise AssertionError(
    "expected parse error for %r (kwargs=%r)" % (input_string[:60], kwargs))


# ─── Strict mode (default) unchanged ───────────────────────────────

def test_default_is_strict():
  _expect_parse_error("_a 1\n")

def test_strict_true_rejects_leading_pair():
  _expect_parse_error("_a 1\ndata_foo\n_x 1\n", strict=True)


# ─── Non-strict: synthesize global_ block ──────────────────────────

def test_non_strict_leading_pair_then_data_block():
  doc = xcif_ext.parse("_a 1\ndata_foo\n_x 1\n", strict=False)
  assert len(doc) == 2, len(doc)
  assert doc[0].name == "global_", doc[0].name
  assert doc[0].find_value("_a") == "1"
  assert doc[1].name == "foo"
  assert doc[1].find_value("_x") == "1"

def test_non_strict_only_leading_pairs():
  doc = xcif_ext.parse("_a 1\n_b 2\n", strict=False)
  assert len(doc) == 1
  assert doc[0].name == "global_"
  assert doc[0].find_value("_a") == "1"
  assert doc[0].find_value("_b") == "2"

def test_non_strict_monomer_library_header_pattern():
  src = (
    "_lib_update       15/04/05\n"
    "# ------------------------------------------------\n"
    "#\n"
    "# ---   LIST OF MONOMERS ---\n"
    "#\n"
    "data_comp_list\n"
    "loop_\n"
    "_chem_comp.id\n"
    "ALA\n"
    "GLY\n"
  )
  doc = xcif_ext.parse(src, strict=False)
  assert len(doc) == 2
  assert doc[0].name == "global_"
  assert doc[0].find_value("_lib_update") == "15/04/05"
  assert doc[1].name == "comp_list"
  lp = doc[1].find_loop("_chem_comp.id")
  assert lp is not None
  assert lp.length == 2

def test_non_strict_leading_loop():
  src = "loop_\n_a\n_b\n1 2\n3 4\n"
  doc = xcif_ext.parse(src, strict=False)
  assert len(doc) == 1
  assert doc[0].name == "global_"
  lp = doc[0].find_loop("_a")
  assert lp is not None
  assert lp.length == 2
  assert lp.width == 2

def test_non_strict_empty_input():
  doc = xcif_ext.parse("", strict=False)
  assert len(doc) == 0

def test_non_strict_only_data_block_unchanged():
  doc = xcif_ext.parse("data_foo\n_x 1\n", strict=False)
  assert len(doc) == 1
  assert doc[0].name == "foo"

def test_non_strict_find_block_by_global_name():
  doc = xcif_ext.parse("_a 1\ndata_foo\n", strict=False)
  g = doc.find_block("global_")
  assert g is not None
  assert g.name == "global_"


# ─── Explicit global_ block header (CIF 1.1 reserved) ─────────────

def test_explicit_global_header_strict():
  doc = xcif_ext.parse("global_\n_a 1\n_b 2\n")
  assert len(doc) == 1
  assert doc[0].name == "global_"
  assert doc[0].find_value("_a") == "1"
  assert doc[0].find_value("_b") == "2"

def test_explicit_global_followed_by_data_block():
  doc = xcif_ext.parse("global_\n_a 1\ndata_foo\n_x 1\n")
  assert len(doc) == 2
  assert doc[0].name == "global_"
  assert doc[0].find_value("_a") == "1"
  assert doc[1].name == "foo"

def test_explicit_global_case_insensitive():
  doc = xcif_ext.parse("GLOBAL_\n_a 1\n")
  assert len(doc) == 1
  g = doc.find_block("global_")
  assert g is not None

def test_monomer_library_header_with_explicit_global():
  src = (
    "global_\n"
    "_lib_name         mon_lib\n"
    "_lib_version      4.11\n"
    "_lib_update       15/04/05\n"
    "# -----------------------------------------------\n"
    "data_comp_list\n"
    "loop_\n"
    "_chem_comp.id\n"
    "ALA\n"
  )
  doc = xcif_ext.parse(src)  # strict=True
  assert len(doc) == 2
  assert doc[0].name == "global_"
  assert doc[0].find_value("_lib_version") == "4.11"
  assert doc[1].name == "comp_list"


def run():
  test_default_is_strict()
  test_strict_true_rejects_leading_pair()
  test_non_strict_leading_pair_then_data_block()
  test_non_strict_only_leading_pairs()
  test_non_strict_monomer_library_header_pattern()
  test_non_strict_leading_loop()
  test_non_strict_empty_input()
  test_non_strict_only_data_block_unchanged()
  test_non_strict_find_block_by_global_name()
  test_explicit_global_header_strict()
  test_explicit_global_followed_by_data_block()
  test_explicit_global_case_insensitive()
  test_monomer_library_header_with_explicit_global()
  print("OK")


if __name__ == "__main__":
  run()
