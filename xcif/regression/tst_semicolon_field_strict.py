"""CIF 1.1 §2.2.7.4: the closing delimiter of a semicolon text field
must be a `;` at column 1 of a new line. An unterminated field (no
such closer before EOF) is a syntax error.

Previously xcif returned a partial TOKEN_VALUE for unterminated fields
and the parse succeeded silently — masking real malformations.
"""
from __future__ import absolute_import, division, print_function

import xcif_ext


def _expect_parse_error(input_string, **kwargs):
  try:
    xcif_ext.parse(input_string, **kwargs)
  except (RuntimeError, ValueError):
    return
  raise AssertionError(
    "expected parse error for %r" % input_string[:60])


# ─── Unterminated fields are errors ───────────────────────────────

def test_semicolon_field_eof_before_close():
  _expect_parse_error("data_t\n_a\n;line one\nline two\n")

def test_semicolon_field_with_mid_line_semicolon():
  # Mirrors tst_lex_parse_build.py::bad_semicolon_text_field
  src = (
    "data_sucrose\n"
    "_a 1\n"
    "_exptl_absorpt_process_details\n"
    ";\n"
    "Final HKLF 4 output contains 64446 reflections, Rint = 0.0650\n"
    " (47528 with I > 3sig(I), Rint = 0.0624);\n"
  )
  _expect_parse_error(src)

def test_semicolon_field_eof_at_open_line():
  _expect_parse_error("data_t\n_a\n;no closing")


# ─── Well-formed semicolon fields continue to parse ───────────────

def test_semicolon_field_properly_closed():
  doc = xcif_ext.parse("data_t\n_a\n;hello\n;\n")
  assert len(doc) == 1
  assert doc[0].find_value("_a") == "hello\n"

def test_semicolon_field_multiline_closed():
  doc = xcif_ext.parse(
    "data_t\n"
    "_a\n"
    ";line one\n"
    "line two\n"
    "line three\n"
    ";\n")
  assert len(doc) == 1

def test_semicolon_field_empty():
  doc = xcif_ext.parse("data_t\n_a\n;\n;\n")
  assert len(doc) == 1


def run():
  test_semicolon_field_eof_before_close()
  test_semicolon_field_with_mid_line_semicolon()
  test_semicolon_field_eof_at_open_line()
  test_semicolon_field_properly_closed()
  test_semicolon_field_multiline_closed()
  test_semicolon_field_empty()
  print("OK")


if __name__ == "__main__":
  run()
