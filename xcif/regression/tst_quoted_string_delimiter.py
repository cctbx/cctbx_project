"""CIF 1.1 syntax spec, paragraph 15
(https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax): a
quote-delimited string may contain instances of the delimiter
provided they are not followed by whitespace. The closing `'` (or
`"`) is only recognised when followed by whitespace or EOL; an
embedded delimiter followed by a non-whitespace char is part of
the string content.

Motivating case — cctbx monomer library, mon_lib_list.cif line 1688:
  'N-METHYL-PYRIDOXAL-5'-PHOSPHATE     '
The inner `'` is followed by `-`, so it must not terminate the string.
"""
from __future__ import absolute_import, division, print_function

import xcif_ext


# ─── Existing behaviour still works ───────────────────────────────

def test_plain_single_quoted():
  doc = xcif_ext.parse("data_t\n_a 'hello'\n")
  assert doc[0].find_value("_a") == "hello"

def test_plain_double_quoted():
  doc = xcif_ext.parse('data_t\n_a "hello"\n')
  assert doc[0].find_value("_a") == "hello"


# ─── Embedded delimiter (CIF 1.1 spec requirement) ────────────────

def test_apostrophe_followed_by_dash_is_embedded():
  doc = xcif_ext.parse("data_t\n_a 'N-METHYL-PYRIDOXAL-5'-PHOSPHATE'\n")
  assert doc[0].find_value("_a") == "N-METHYL-PYRIDOXAL-5'-PHOSPHATE"

def test_apostrophe_followed_by_letter_is_embedded():
  doc = xcif_ext.parse("data_t\n_a 'it's fine'\n")
  assert doc[0].find_value("_a") == "it's fine"

def test_double_quote_followed_by_letter_is_embedded():
  # No space before `once` — else the `" ` sequence would close the
  # string per CIF 1.1 paragraph 15.
  doc = xcif_ext.parse('data_t\n_a "say "hi"once"\n')
  assert doc[0].find_value("_a") == 'say "hi"once'

def test_multiple_embedded_apostrophes():
  doc = xcif_ext.parse("data_t\n_a 'don't can't won't'\n")
  assert doc[0].find_value("_a") == "don't can't won't"


# ─── Closing rules — whitespace / EOL / tab ───────────────────────

def test_apostrophe_at_end_of_input_closes():
  doc = xcif_ext.parse("data_t\n_a 'hello'")
  assert doc[0].find_value("_a") == "hello"

def test_apostrophe_followed_by_tab_closes():
  # Tab after the closing `'` must be recognised as whitespace so the
  # string closes and the next tag is recognised.
  doc = xcif_ext.parse("data_t\n_a 'hello'\t_b 3\n")
  assert doc[0].find_value("_a") == "hello"
  assert doc[0].find_value("_b") == "3"


# ─── Monomer library real-world case ──────────────────────────────

def test_monomer_library_apostrophe_row():
  src = (
    "data_t\n"
    "loop_\n"
    "_chem_comp.id\n"
    "_chem_comp.three_letter_code\n"
    "_chem_comp.name\n"
    "_chem_comp.group\n"
    "_chem_comp.number_atoms_all\n"
    "_chem_comp.number_atoms_nh\n"
    "_chem_comp.desc_level\n"
    "MPL      .   'N-METHYL-PYRIDOXAL-5'-PHOSPHATE     ' "
    "non-polymer  28  17 M\n"
  )
  doc = xcif_ext.parse(src)
  lp = doc[0].find_loop("_chem_comp.id")
  assert lp is not None
  assert lp.length == 1
  assert lp.width == 7
  assert lp.value(0, 0) == "MPL"
  assert lp.value(0, 2) == "N-METHYL-PYRIDOXAL-5'-PHOSPHATE     "
  assert lp.value(0, 3) == "non-polymer"
  assert lp.value(0, 6) == "M"


def run():
  test_plain_single_quoted()
  test_plain_double_quoted()
  test_apostrophe_followed_by_dash_is_embedded()
  test_apostrophe_followed_by_letter_is_embedded()
  test_double_quote_followed_by_letter_is_embedded()
  test_multiple_embedded_apostrophes()
  test_apostrophe_at_end_of_input_closes()
  test_apostrophe_followed_by_tab_closes()
  test_monomer_library_apostrophe_row()
  print("OK")


if __name__ == "__main__":
  run()
