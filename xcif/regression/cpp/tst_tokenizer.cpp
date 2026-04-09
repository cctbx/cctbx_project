// cctbx_project/xcif/regression/cpp/tst_tokenizer.cpp
//
// TDD tokenizer tests — Categories A (basic tokens) and B (edge cases).
// These will NOT compile until xcif/tokenizer.h is implemented (Step 4).
// See test_plan_2.txt §3.1 and §3.2.

#include "test_utils.h"
#include "xcif/tokenizer.h"
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Tokenize all tokens from a NUL-terminated CIF string.
static std::vector<xcif::Token> tokenize_all(const char* input) {
  xcif::Tokenizer tok(input, std::strlen(input), "<test>");
  std::vector<xcif::Token> out;
  for (;;) {
    xcif::Token t = tok.next();
    out.push_back(t);
    if (t.type == xcif::TOKEN_EOF) break;
  }
  return out;
}

// Return the Nth token (0-indexed).  Returns TOKEN_EOF if input is shorter.
static xcif::Token nth(const char* input, int n) {
  xcif::Tokenizer tok(input, std::strlen(input), "<test>");
  xcif::Token t;
  for (int i = 0; i <= n; ++i) t = tok.next();
  return t;
}

static xcif::Token first(const char* s)  { return nth(s, 0); }
static xcif::Token second(const char* s) { return nth(s, 1); }
static xcif::Token third(const char* s)  { return nth(s, 2); }

// ---------------------------------------------------------------------------
// Category A: Basic Token Types (§3.1)
// ---------------------------------------------------------------------------

static void test_eof_on_empty() {
  CHECK_EQ(first("").type, xcif::TOKEN_EOF);
}

static void test_eof_on_whitespace_only() {
  CHECK_EQ(first("   \t  \n  ").type, xcif::TOKEN_EOF);
}

static void test_tag_token() {
  xcif::Token t = first("_cell.length_a");
  CHECK_EQ(t.type, xcif::TOKEN_TAG);
  CHECK_EQ(t.as_str(), std::string("_cell.length_a"));
}

static void test_unquoted_value() {
  xcif::Token t = second("_a 5.432");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("5.432"));
}

static void test_single_quoted_string() {
  xcif::Token t = second("_a 'hello world'");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("hello world"));
}

static void test_double_quoted_string() {
  xcif::Token t = second("_a \"hello world\"");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("hello world"));
}

static void test_semicolon_text_field() {
  // Value = everything between the newline after ';' and the line-initial ';'
  const char* input = "_a\n;line one\nline two\n;\n";
  xcif::Token t = second(input);
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("line one\nline two\n"));
}

static void test_cif2_triple_double_quoted() {
  xcif::Token t = second("_a \"\"\"triple quoted\"\"\"");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("triple quoted"));
}

static void test_cif2_triple_single_quoted() {
  xcif::Token t = second("_a '''triple quoted'''");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("triple quoted"));
}

static void test_block_header() {
  xcif::Token t = first("data_1YJP");
  CHECK_EQ(t.type, xcif::TOKEN_BLOCK_HEADER);
  CHECK_EQ(t.as_str(), std::string("data_1YJP"));
}

static void test_loop_keyword() {
  CHECK_EQ(first("loop_").type, xcif::TOKEN_LOOP);
}

static void test_save_frame_header() {
  xcif::Token t = first("save_FOO");
  CHECK_EQ(t.type, xcif::TOKEN_SAVE_HEADER);
  CHECK_EQ(t.as_str(), std::string("save_FOO"));
}

static void test_save_frame_end() {
  // save_ with nothing after it (or only whitespace) ends a save frame
  CHECK_EQ(first("save_").type, xcif::TOKEN_SAVE_END);
  CHECK_EQ(first("save_ ").type, xcif::TOKEN_SAVE_END);
}

static void test_integer_value() {
  xcif::Token t = second("_a 42");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("42"));
}

static void test_float_value() {
  xcif::Token t = second("_a 3.14159");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("3.14159"));
}

static void test_scientific_notation() {
  xcif::Token t = second("_a 1.23e-10");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("1.23e-10"));
}

static void test_standard_uncertainty() {
  // SU notation: the whole "1.542(3)" is one value token — caller interprets it
  xcif::Token t = second("_a 1.542(3)");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("1.542(3)"));
}

static void test_unknown_value_dot() {
  xcif::Token t = second("_a .");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("."));
}

static void test_inapplicable_value_question() {
  xcif::Token t = second("_a ?");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("?"));
}

static void test_comment_is_skipped() {
  // Comment tokens must be invisible; only TAG VALUE TAG VALUE EOF appear
  std::vector<xcif::Token> toks =
    tokenize_all("_a 1 # ignored comment\n_b 2");
  CHECK_EQ((int)toks.size(), 5);
  CHECK_EQ(toks[0].type, xcif::TOKEN_TAG);
  CHECK_EQ(toks[1].type, xcif::TOKEN_VALUE);
  CHECK_EQ(toks[2].type, xcif::TOKEN_TAG);
  CHECK_EQ(toks[3].type, xcif::TOKEN_VALUE);
  CHECK_EQ(toks[4].type, xcif::TOKEN_EOF);
}

static void test_multiple_blocks() {
  std::vector<xcif::Token> toks =
    tokenize_all("data_A _x 1\ndata_B _y 2");
  // BLOCK TAG VALUE BLOCK TAG VALUE EOF = 7
  CHECK_EQ((int)toks.size(), 7);
  CHECK_EQ(toks[0].type, xcif::TOKEN_BLOCK_HEADER);
  CHECK_EQ(toks[3].type, xcif::TOKEN_BLOCK_HEADER);
}

// ---------------------------------------------------------------------------
// Category B: Edge Cases (§3.2)
// ---------------------------------------------------------------------------

static void test_loop_something_is_value() {
  // "loop_something" — has text after loop_, so NOT the loop_ keyword
  xcif::Token t = second("_a loop_something");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("loop_something"));
}

static void test_data_underscore_as_substring_is_value() {
  // "somedata_block" contains data_ but doesn't start with it → plain value
  xcif::Token t = second("_a somedata_block");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("somedata_block"));
}

static void test_empty_single_quoted_string() {
  xcif::Token t = second("_a ''");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string(""));
}

static void test_empty_double_quoted_string() {
  xcif::Token t = second("_a \"\"");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string(""));
}

static void test_double_quote_inside_single_quoted() {
  xcif::Token t = second("_a 'value with \"quotes\" inside'");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("value with \"quotes\" inside"));
}

static void test_single_quote_inside_double_quoted() {
  xcif::Token t = second("_a \"it's fine\"");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("it's fine"));
}

static void test_crlf_line_endings() {
  // \r\n must be treated as a single newline for line counting
  std::vector<xcif::Token> toks = tokenize_all("_a 1\r\n_b 2");
  CHECK_EQ((int)toks.size(), 5); // TAG VALUE TAG VALUE EOF
  CHECK_EQ(toks[2].as_str(), std::string("_b"));
}

static void test_cr_only_line_endings() {
  std::vector<xcif::Token> toks = tokenize_all("_a 1\r_b 2");
  CHECK_EQ((int)toks.size(), 5);
  CHECK_EQ(toks[2].as_str(), std::string("_b"));
}

static void test_tab_as_whitespace() {
  xcif::Token t = second("_a\t42");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("42"));
}

static void test_utf8_bom_skipped() {
  // UTF-8 BOM (\xEF\xBB\xBF) at file start must be silently consumed
  const char* input = "\xEF\xBB\xBF" "data_test";
  xcif::Token t = first(input);
  CHECK_EQ(t.type, xcif::TOKEN_BLOCK_HEADER);
  CHECK_EQ(t.as_str(), std::string("data_test"));
}

static void test_line_column_of_first_tag() {
  xcif::Token t = first("_cell.a 5.0");
  CHECK_EQ(t.line, 1);
  CHECK_EQ(t.col,  1);
}

static void test_line_column_of_value_after_tag() {
  xcif::Token t = second("_cell.a 5.0");
  CHECK_EQ(t.line, 1);
  CHECK_EQ(t.col,  9); // "_cell.a " = 8 chars, value starts at col 9
}

static void test_line_number_increments_on_newline() {
  xcif::Tokenizer tok("_a 1\n_b 2", 9, "<test>");
  tok.next(); // _a
  tok.next(); // 1
  xcif::Token b = tok.next(); // _b  (line 2)
  CHECK_EQ(b.type, xcif::TOKEN_TAG);
  CHECK_EQ(b.line, 2);
  CHECK_EQ(b.col,  1);
}

static void test_semicolon_not_at_line_start_is_value_char() {
  // A semicolon mid-line is not a text-block delimiter
  const char* input = "_a\n;value; more\n;\n";
  xcif::Token t = second(input);
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("value; more\n"));
}

static void test_no_trailing_newline() {
  // File without trailing newline must still tokenize the last token
  xcif::Token t = first("data_test");
  CHECK_EQ(t.type, xcif::TOKEN_BLOCK_HEADER);
  CHECK_EQ(t.as_str(), std::string("data_test"));
}

static void test_multiple_blank_lines_ignored() {
  xcif::Token t = second("_a\n\n\n\n42");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("42"));
}

static void test_consecutive_tags_and_values() {
  std::vector<xcif::Token> toks =
    tokenize_all("_a 1 _b 2 _c 3");
  // TAG VALUE TAG VALUE TAG VALUE EOF = 7
  CHECK_EQ((int)toks.size(), 7);
  for (int i = 0; i < 6; i += 2) CHECK_EQ(toks[i].type, xcif::TOKEN_TAG);
  for (int i = 1; i < 6; i += 2) CHECK_EQ(toks[i].type, xcif::TOKEN_VALUE);
}

static void test_loop_followed_by_tags_and_values() {
  std::vector<xcif::Token> toks =
    tokenize_all("loop_ _a _b 1 2 3 4");
  CHECK_EQ(toks[0].type, xcif::TOKEN_LOOP);
  CHECK_EQ(toks[1].type, xcif::TOKEN_TAG);
  CHECK_EQ(toks[2].type, xcif::TOKEN_TAG);
  CHECK_EQ(toks[3].type, xcif::TOKEN_VALUE);
}

static void test_unterminated_single_quote_returns_rest() {
  // Unterminated single-quoted string: returns everything after the opening
  // quote up to EOF as TOKEN_VALUE (graceful degradation, not an error).
  xcif::Token t = second("_a 'hello world");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("hello world"));
}

static void test_unterminated_double_quote_returns_rest() {
  xcif::Token t = second("_a \"hello world");
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("hello world"));
}

static void test_null_byte_truncates_unquoted_token() {
  // Null byte is not an ordinary character, so read_unquoted stops there.
  const char input[] = {'_', 'a', ' ', 'v', 'a', 'l', '\0', 'u', 'e'};
  xcif::Tokenizer tok(input, sizeof(input), "<test>");
  tok.next(); // _a
  xcif::Token t = tok.next();
  CHECK_EQ(t.type, xcif::TOKEN_VALUE);
  CHECK_EQ(t.as_str(), std::string("val"));
}

// ---------------------------------------------------------------------------
// Test Registry
// ---------------------------------------------------------------------------

static void run_all_tests() {
  // Category A
  // 3.1 Category A: Basic Token Types

  // Pure C++ tests (standalone executables, no Python). Each test provides a minimal CIF byte buffer and
  // asserts the tokenizer emits the expected token sequence.

  // Token types to cover: unquoted values, single-quoted strings, double-quoted strings, semicolon text fields,
  //  CIF2 triple-quoted strings ("""...""" and '''...'''), data_xxxx block headers, save_xxxx and save_ frame
  // headers/terminators, loop_, tag names (_tag.name), integers, floats, scientific notation, values with
  // standard uncertainties (1.542(3)), unknown value (.), inapplicable value (?), and comments (# ...).

  test_eof_on_empty();
  test_eof_on_whitespace_only();
  test_tag_token();
  test_unquoted_value();
  test_single_quoted_string();
  test_double_quoted_string();
  test_semicolon_text_field();
  test_cif2_triple_double_quoted();
  test_cif2_triple_single_quoted();
  test_block_header();
  test_loop_keyword();
  test_save_frame_header();
  test_save_frame_end();
  test_integer_value();
  test_float_value();
  test_scientific_notation();
  test_standard_uncertainty();
  test_unknown_value_dot();
  test_inapplicable_value_question();
  test_comment_is_skipped();
  test_multiple_blocks();
  // Category B
  // 3.2 Category B: Tokenizer Edge Cases

  // Pure C++ tests. Targets ambiguous or unusual input a naive tokenizer might mishandle.

  // Must cover: values that resemble keywords but are not (loop_something as a data value, data_ as a
  // substring), empty quoted strings ('', ""), semicolon text fields containing lines beginning with semicolons
  //  (leading backslash or space per spec), semicolon text fields with trailing whitespace on the delimiter
  // line, nested quotes ('it''s' in CIF2, "value with 'quotes'"), very long unquoted values (10,000+
  // characters), tag names with unusual but legal characters, Unicode content in CIF2, mixed line endings (\n,
  // \r\n, \r), files with and without trailing newline, byte order marks (BOM) at file start, tab characters as
  //  whitespace, consecutive whitespace and multiple blank lines.

  test_loop_something_is_value();
  test_data_underscore_as_substring_is_value();
  test_empty_single_quoted_string();
  test_empty_double_quoted_string();
  test_double_quote_inside_single_quoted();
  test_single_quote_inside_double_quoted();
  test_crlf_line_endings();
  test_cr_only_line_endings();
  test_tab_as_whitespace();
  test_utf8_bom_skipped();
  test_line_column_of_first_tag();
  test_line_column_of_value_after_tag();
  test_line_number_increments_on_newline();
  test_semicolon_not_at_line_start_is_value_char();
  test_no_trailing_newline();
  test_multiple_blank_lines_ignored();
  test_consecutive_tags_and_values();
  test_loop_followed_by_tags_and_values();
  test_unterminated_single_quote_returns_rest();
  test_unterminated_double_quote_returns_rest();
  test_null_byte_truncates_unquoted_token();
}

XCIF_TEST_MAIN()
