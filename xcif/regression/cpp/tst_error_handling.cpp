// cctbx_project/xcif/regression/cpp/tst_error_handling.cpp
#include <xcif/data_model.h>
#include "test_utils.h"
#include <string>
#include <stdexcept>

using namespace xcif;

// Helper: check that parsing throws CifError containing expected substring
static bool parse_throws(const char* input, const char* expected_substr) {
  try {
    parse(input);
    return false;
  } catch (const CifError& e) {
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

// Helper: check that parse throws CifError at expected line
static bool parse_throws_at_line(const char* input, int expected_line) {
  try {
    parse(input);
    return false;
  } catch (const CifError& e) {
    if (e.line() != expected_line) {
      fprintf(stderr, "  Expected error at line %d, got line %d\n",
              expected_line, e.line());
      return false;
    }
    return true;
  } catch (const std::exception&) {
    return false;
  }
}

// ─── §1  Structural errors ─────────────────────────────────────────

void test_tag_outside_block() {
  CHECK(parse_throws("_tag value\n", "outside"));
}

void test_value_outside_block() {
  CHECK(parse_throws("some_value\n", "outside"));
}

void test_loop_outside_block() {
  CHECK(parse_throws("loop_\n_tag\n1\n", "outside"));
}

void test_tag_missing_value() {
  // Tag followed immediately by another tag
  CHECK(parse_throws("data_t\n_a\n_b val\n", "value"));
}

void test_tag_missing_value_at_eof() {
  CHECK(parse_throws("data_t\n_a\n", "value"));
}

void test_loop_no_tags() {
  // loop_ followed by a value instead of a tag
  CHECK(parse_throws("data_t\nloop_\nvalue\n", "tag"));
}

void test_loop_values_not_multiple_of_tags() {
  // 2 tags but 3 values
  CHECK(parse_throws(
    "data_t\nloop_\n_a\n_b\n1 2 3\n",
    "multiple"
  ));
}

// ─── §2  Duplicate detection ───────────────────────────────────────

void test_duplicate_tag_in_block() {
  CHECK(parse_throws(
    "data_t\n_x 1\n_x 2\n",
    "duplicate"
  ));
}

void test_duplicate_tag_case_insensitive() {
  CHECK(parse_throws(
    "data_t\n_Tag 1\n_TAG 2\n",
    "duplicate"
  ));
}

void test_duplicate_tag_in_loop() {
  CHECK(parse_throws(
    "data_t\nloop_\n_a\n_a\n1 2\n",
    "duplicate"
  ));
}

// ─── §3  Error location tracking ───────────────────────────────────

void test_error_reports_line_number() {
  CHECK(parse_throws_at_line(
    "data_t\n_ok 1\n_bad\n",  // error on line 3 (missing value)
    3
  ));
}

void test_error_reports_source_name() {
  try {
    parse("_orphan val\n", "<test>");
    CHECK(false);  // should have thrown
  } catch (const CifError& e) {
    std::string msg = e.what();
    CHECK(msg.find("<test>") != std::string::npos);
  }
}

// ─── §4  Save frame errors ─────────────────────────────────────────

void test_save_end_without_save_start() {
  CHECK(parse_throws(
    "data_t\nsave_\n",
    "save"
  ));
}

void test_unclosed_save_frame() {
  CHECK(parse_throws(
    "data_t\nsave_frame\n_x 1\n",
    "save"
  ));
}

// ─── §5  Graceful handling ──────────────────────────────────────────

void test_empty_input_no_error() {
  Document doc = parse("");
  CHECK_EQ(doc.size(), 0u);
}

void test_comment_only_no_error() {
  Document doc = parse("# just a comment\n");
  CHECK_EQ(doc.size(), 0u);
}

void test_whitespace_only_no_error() {
  Document doc = parse("   \n  \n  \n");
  CHECK_EQ(doc.size(), 0u);
}

void test_unterminated_quoted_string_graceful() {
  // Missing closing quote: tokenizer returns partial content and parse succeeds.
  Document doc = parse("data_t\n_a 'hello");
  CHECK_EQ(doc.size(), 1u);
  CHECK(doc[0].has_tag(string_view("_a")));
}

void test_unclosed_semicolon_field_raises() {
  // CIF 1.1 §2.2.7.4: semicolon text field must be closed by `;` at
  // column 1 of a new line. EOF before that is a syntax error.
  // (Was previously "graceful"; tightened to match spec + ucif.)
  CHECK(parse_throws("data_t\n_a\n;line one\nline two\n", "semicolon"));
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  // §1 Structural errors
  test_tag_outside_block();
  test_value_outside_block();
  test_loop_outside_block();
  test_tag_missing_value();
  test_tag_missing_value_at_eof();
  test_loop_no_tags();
  test_loop_values_not_multiple_of_tags();

  // §2 Duplicate detection
  test_duplicate_tag_in_block();
  test_duplicate_tag_case_insensitive();
  test_duplicate_tag_in_loop();

  // §3 Error location
  test_error_reports_line_number();
  test_error_reports_source_name();

  // §4 Save frame errors
  test_save_end_without_save_start();
  test_unclosed_save_frame();

  // §5 Graceful handling
  test_empty_input_no_error();
  test_comment_only_no_error();
  test_whitespace_only_no_error();
  test_unterminated_quoted_string_graceful();
  test_unclosed_semicolon_field_raises();
}

XCIF_TEST_MAIN()
