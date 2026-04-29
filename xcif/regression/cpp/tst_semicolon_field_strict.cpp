// cctbx_project/xcif/regression/cpp/tst_semicolon_field_strict.cpp
//
// CIF 1.1 syntax spec, paragraph 17
// (https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax):
// a semicolon text field is delimited by `;` appearing as the first
// character of a line. An unterminated field (EOF before the closing
// `;` at column 1, or only a mid-line `;`) is a syntax error.
//
// xcif previously returned a partial TOKEN_VALUE and let the parse
// succeed, silently accepting content that has no closing delimiter.
// That hid real malformations in user files. This test fixture
// asserts the parser now raises CifError for each malformed shape.

#include <xcif/data_model.h>
#include "test_utils.h"
#include <string>

using namespace xcif;

static bool parse_throws(const char* input, const char* expected_substr) {
  try {
    parse(input);
    return false;
  } catch (const CifError& e) {
    if (expected_substr == nullptr) return true;
    std::string msg = e.what();
    if (msg.find(expected_substr) == std::string::npos) {
      fprintf(stderr, "  Expected substring '%s' in error:\n  '%s'\n",
              expected_substr, msg.c_str());
      return false;
    }
    return true;
  } catch (const std::exception& e) {
    fprintf(stderr, "  Wrong exception type: %s\n", e.what());
    return false;
  }
}

// ─── Unterminated fields are errors ───────────────────────────────

void test_semicolon_field_eof_before_close() {
  // Opens at col 1 but never closes.
  CHECK(parse_throws("data_t\n_a\n;line one\nline two\n",
                     "semicolon"));
}

void test_semicolon_field_with_mid_line_semicolon() {
  // The trailing `;` is at end-of-line, not column 1 — must not close.
  // This is the shape used by tst_lex_parse_build.py::bad_semicolon_text_field.
  const char* src =
    "data_sucrose\n"
    "_a 1\n"
    "_exptl_absorpt_process_details\n"
    ";\n"
    "Final HKLF 4 output contains 64446 reflections, Rint = 0.0650\n"
    " (47528 with I > 3sig(I), Rint = 0.0624);\n";
  CHECK(parse_throws(src, "semicolon"));
}

void test_semicolon_field_eof_at_open_line() {
  // File ends on the opening `;` line itself.
  CHECK(parse_throws("data_t\n_a\n;no closing", "semicolon"));
}

// ─── Well-formed semicolon fields continue to parse ───────────────

void test_semicolon_field_properly_closed_two_lines() {
  Document doc = parse("data_t\n_a\n;hello\n;\n");
  CHECK_EQ(doc.size(), 1u);
  CHECK(doc[0].has_tag(string_view("_a")));
  // Value is the body between the opening `;\n` and the closing `\n;`.
  // Content begins immediately after the opening `;`, so includes an
  // initial newline.
  string_view v = doc[0].find_value(string_view("_a"));
  CHECK_EQ(std::string(v.data(), v.size()), std::string("hello\n"));
}

void test_semicolon_field_multiline_closed() {
  Document doc = parse(
    "data_t\n"
    "_a\n"
    ";line one\n"
    "line two\n"
    "line three\n"
    ";\n");
  CHECK_EQ(doc.size(), 1u);
  CHECK(doc[0].has_tag(string_view("_a")));
}

void test_semicolon_field_empty() {
  // Immediate close on the next line.
  Document doc = parse("data_t\n_a\n;\n;\n");
  CHECK_EQ(doc.size(), 1u);
  CHECK(doc[0].has_tag(string_view("_a")));
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  test_semicolon_field_eof_before_close();
  test_semicolon_field_with_mid_line_semicolon();
  test_semicolon_field_eof_at_open_line();
  test_semicolon_field_properly_closed_two_lines();
  test_semicolon_field_multiline_closed();
  test_semicolon_field_empty();
}

XCIF_TEST_MAIN()
