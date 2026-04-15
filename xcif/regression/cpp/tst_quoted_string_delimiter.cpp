// cctbx_project/xcif/regression/cpp/tst_quoted_string_delimiter.cpp
//
// CIF 1.1 syntax spec, paragraph 15
// (https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax):
// a quote-delimited character string may contain instances of the
// delimiter (' for single, " for double) provided they are not
// followed by whitespace. The closing delimiter is only recognised
// as such when it is followed by whitespace or end-of-input.
//
// Motivating case — cctbx monomer library, mon_lib_list.cif line 1688:
//   'N-METHYL-PYRIDOXAL-5'-PHOSPHATE     '
// The inner `'` is followed by `-`, so the string runs to the outer
// closing `'` which is followed by a space.

#include <xcif/data_model.h>
#include <xcif/tokenizer.h>
#include "test_utils.h"
#include <cstring>
#include <string>

using namespace xcif;

static std::string sv_to_string(const string_view& sv) {
  return std::string(sv.data(), sv.size());
}

// Tokenize once and return the Nth value token (0-indexed).
static std::string nth_value(const char* input, std::size_t n) {
  Tokenizer tok(input, std::strlen(input), "<test>");
  std::size_t seen = 0;
  while (true) {
    Token t = tok.next();
    if (t.type == TOKEN_EOF) break;
    if (t.type == TOKEN_VALUE) {
      if (seen == n) return std::string(t.ptr, t.len);
      ++seen;
    }
  }
  return std::string("<not found>");
}

// ─── Existing behaviour still works ────────────────────────────────

void test_plain_single_quoted() {
  CHECK_EQ(nth_value("_a 'hello'\n", 0), std::string("hello"));
}

void test_plain_double_quoted() {
  CHECK_EQ(nth_value("_a \"hello\"\n", 0), std::string("hello"));
}

void test_empty_single_quoted() {
  CHECK_EQ(nth_value("_a '' \n", 0), std::string(""));
}

// ─── Embedded delimiter (CIF 1.1 spec requirement) ────────────────

void test_apostrophe_followed_by_dash_is_embedded() {
  // `5'-PHOSPHATE` — apostrophe followed by `-`, not whitespace.
  CHECK_EQ(nth_value("_a 'N-METHYL-PYRIDOXAL-5'-PHOSPHATE'\n", 0),
           std::string("N-METHYL-PYRIDOXAL-5'-PHOSPHATE"));
}

void test_apostrophe_followed_by_letter_is_embedded() {
  CHECK_EQ(nth_value("_a 'it's fine'\n", 0),
           std::string("it's fine"));
}

void test_double_quote_followed_by_letter_is_embedded() {
  // Same rule applies to double-quoted strings (paragraph 15).
  // Note the absence of a space between the inner `"` and `once` —
  // a space there would close the string.
  CHECK_EQ(nth_value("_a \"say \"hi\"once\"\n", 0),
           std::string("say \"hi\"once"));
}

void test_multiple_embedded_apostrophes() {
  CHECK_EQ(nth_value("_a 'don't can't won't' \n", 0),
           std::string("don't can't won't"));
}

// ─── Closing rules — whitespace / EOL / tab ───────────────────────

void test_apostrophe_at_end_of_line_closes() {
  // Closing `'` followed by newline is a terminator.
  CHECK_EQ(nth_value("_a 'hello'\n_b 1\n", 0), std::string("hello"));
}

void test_apostrophe_followed_by_tab_closes() {
  CHECK_EQ(nth_value("_a 'hello'\t2\n", 0), std::string("hello"));
}

void test_apostrophe_at_end_of_input_closes() {
  // No character after the closing `'`.
  CHECK_EQ(nth_value("_a 'hello'", 0), std::string("hello"));
}

// ─── Monomer library real-world case ──────────────────────────────

void test_monomer_library_apostrophe_row() {
  // Parse a compact loop row like the one at mon_lib_list.cif:1688.
  const char* src =
    "data_t\n"
    "loop_\n"
    "_chem_comp.id\n"
    "_chem_comp.three_letter_code\n"
    "_chem_comp.name\n"
    "_chem_comp.group\n"
    "_chem_comp.number_atoms_all\n"
    "_chem_comp.number_atoms_nh\n"
    "_chem_comp.desc_level\n"
    "MPL      .   'N-METHYL-PYRIDOXAL-5'-PHOSPHATE     ' non-polymer  28  17 M\n";
  Document doc = parse(src);
  CHECK_EQ(doc.size(), 1u);
  const Loop* lp = doc[0].find_loop(string_view("_chem_comp.id"));
  CHECK(lp != 0);
  if (!lp) return;
  CHECK_EQ(lp->length(), 1u);
  CHECK_EQ(lp->width(), 7u);
  // Column 2 (0-indexed) is the name with embedded apostrophe.
  CHECK_EQ(sv_to_string(lp->value(0, 2)),
           std::string("N-METHYL-PYRIDOXAL-5'-PHOSPHATE     "));
  // Surrounding columns must be unaffected.
  CHECK_EQ(sv_to_string(lp->value(0, 0)), std::string("MPL"));
  CHECK_EQ(sv_to_string(lp->value(0, 3)), std::string("non-polymer"));
  CHECK_EQ(sv_to_string(lp->value(0, 6)), std::string("M"));
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  test_plain_single_quoted();
  test_plain_double_quoted();
  test_empty_single_quoted();

  test_apostrophe_followed_by_dash_is_embedded();
  test_apostrophe_followed_by_letter_is_embedded();
  test_double_quote_followed_by_letter_is_embedded();
  test_multiple_embedded_apostrophes();

  test_apostrophe_at_end_of_line_closes();
  test_apostrophe_followed_by_tab_closes();
  test_apostrophe_at_end_of_input_closes();

  test_monomer_library_apostrophe_row();
}

XCIF_TEST_MAIN()
