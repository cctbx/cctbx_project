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
  // 2 tags but 3 values. The error message must name the offending
  // loop's first tag so users can grep for it in large files, and
  // include the actual counts so the user understands the shape.
  CHECK(parse_throws(
    "data_t\nloop_\n_alpha\n_beta\n1 2 3\n",
    "_alpha"));
  CHECK(parse_throws(
    "data_t\nloop_\n_alpha\n_beta\n1 2 3\n",
    "3 values"));
  CHECK(parse_throws(
    "data_t\nloop_\n_alpha\n_beta\n1 2 3\n",
    "2 tags"));
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

void test_unterminated_quoted_string_raises() {
  // CIF 1.1 grammar requires a matching closing delimiter for a
  // quoted string; EOF before the close is a syntax error.
  // (See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax.)
  // Was previously "graceful"; tightened to match spec + ucif.
  CHECK(parse_throws("data_t\n_a 'hello", "unterminated"));
  CHECK(parse_throws("data_t\n_a \"hello", "unterminated"));
}

void test_dos_eof_marker_accepted_at_end() {
  // Legacy DOS end-of-file marker (0x1A) commonly appended by SHELX
  // and other tools to `.hkl` / `.cif` files. CIF 1.1 doesn't list
  // it in the character set, but ucif ignores it and real files
  // contain it. xcif must treat it as EOF, not as a stray value
  // that bumps the loop's value count.
  const char src[] =
    "data_t\n"
    "loop_\n_a\n_b\n1 2\n3 4\n"
    "\x1A\n";
  Document doc = parse(src, sizeof(src) - 1);
  CHECK_EQ(doc.size(), 1u);
  const Loop* lp = doc[0].find_loop(string_view("_a"));
  CHECK(lp != 0);
  if (lp) CHECK_EQ(lp->length(), 2u);
}

void test_dos_eof_marker_mid_stream_is_eof() {
  // If 0x1A appears mid-stream (followed by more content), treat it
  // as EOF and ignore everything after — matches DOS-era semantics.
  const char src[] = "data_t\n_a 1\n\x1A\n_b 2\n";
  Document doc = parse(src, sizeof(src) - 1);
  CHECK_EQ(doc.size(), 1u);
  CHECK_EQ(std::string(doc[0].find_value(string_view("_a")).data(),
                       doc[0].find_value(string_view("_a")).size()),
           std::string("1"));
  // _b never appeared — it's past the DOS-EOF
  CHECK(!doc[0].has_tag(string_view("_b")));
}

void test_unclosed_semicolon_field_raises() {
  // CIF 1.1 syntax spec paragraph 17
  // (https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax):
  // a semicolon text field is closed by `;` as the first character of
  // a line. EOF before that is a syntax error.
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
  test_unterminated_quoted_string_raises();
  test_dos_eof_marker_accepted_at_end();
  test_dos_eof_marker_mid_stream_is_eof();
  test_unclosed_semicolon_field_raises();
}

XCIF_TEST_MAIN()
