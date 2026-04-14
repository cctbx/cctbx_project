// cctbx_project/xcif/regression/cpp/tst_non_strict_mode.cpp
//
// Tests for xcif::parse(..., strict=false).
//
// Motivating case: the cctbx monomer library files start with
//   _lib_update       15/04/05
//   # ---   LIST OF MONOMERS ---
//   data_comp_list
//   ...
// i.e. pair items appear before the first data_ block. ucif accepts this
// silently when strict=false; xcif rejected it with a "data outside of a
// data block" error. Non-strict mode synthesizes an implicit "global_"
// block to hold any leading / interleaved content.

#include <xcif/data_model.h>
#include "test_utils.h"
#include <string>

using namespace xcif;

static std::string sv_to_string(const string_view& sv) {
  return std::string(sv.data(), sv.size());
}

// ─── Strict mode unchanged ─────────────────────────────────────────

void test_strict_rejects_leading_pair() {
  try {
    parse("_a 1\ndata_foo\n_x 1\n");
    CHECK(false /* should have thrown */);
  } catch (const CifError&) {
    // expected
  }
}

void test_default_is_strict() {
  // No strict arg — must still error as before.
  try {
    parse("_a 1\n");
    CHECK(false /* should have thrown */);
  } catch (const CifError& e) {
    CHECK(std::string(e.what()).find("outside") != std::string::npos);
  }
}

void test_strict_true_rejects_leading_pair() {
  // Explicit strict=true.
  try {
    parse("_a 1\n", "<test>", /*strict=*/true);
    CHECK(false /* should have thrown */);
  } catch (const CifError&) {
    // expected
  }
}

// ─── Non-strict: synthesize global_ block ─────────────────────────

void test_non_strict_leading_pair_then_data_block() {
  Document doc = parse(
    "_a 1\ndata_foo\n_x 1\n",
    "<test>",
    /*strict=*/false);
  CHECK_EQ(doc.size(), 2u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_a"))),
           std::string("1"));
  CHECK_EQ(sv_to_string(doc[1].name()), std::string("foo"));
  CHECK_EQ(sv_to_string(doc[1].find_value(string_view("_x"))),
           std::string("1"));
}

void test_non_strict_only_leading_pairs() {
  Document doc = parse("_a 1\n_b 2\n", "<test>", false);
  CHECK_EQ(doc.size(), 1u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_a"))),
           std::string("1"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_b"))),
           std::string("2"));
}

void test_non_strict_monomer_library_header_pattern() {
  const char* src =
    "_lib_update       15/04/05\n"
    "# ------------------------------------------------\n"
    "#\n"
    "# ---   LIST OF MONOMERS ---\n"
    "#\n"
    "data_comp_list\n"
    "loop_\n"
    "_chem_comp.id\n"
    "ALA\n"
    "GLY\n";
  Document doc = parse(src, "<monomer>", false);
  CHECK_EQ(doc.size(), 2u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_lib_update"))),
           std::string("15/04/05"));
  CHECK_EQ(sv_to_string(doc[1].name()), std::string("comp_list"));
  const Loop* lp = doc[1].find_loop(string_view("_chem_comp.id"));
  CHECK(lp != 0);
  if (lp) CHECK_EQ(lp->length(), 2u);
}

void test_non_strict_leading_loop() {
  const char* src =
    "loop_\n"
    "_a\n"
    "_b\n"
    "1 2\n"
    "3 4\n";
  Document doc = parse(src, "<test>", false);
  CHECK_EQ(doc.size(), 1u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  const Loop* lp = doc[0].find_loop(string_view("_a"));
  CHECK(lp != 0);
  if (lp) {
    CHECK_EQ(lp->length(), 2u);
    CHECK_EQ(lp->width(), 2u);
  }
}

void test_non_strict_empty_input() {
  Document doc = parse("", "<test>", false);
  CHECK_EQ(doc.size(), 0u);
}

void test_non_strict_only_data_block_unchanged() {
  // If there's no leading content, non-strict should behave like strict.
  Document doc = parse("data_foo\n_x 1\n", "<test>", false);
  CHECK_EQ(doc.size(), 1u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("foo"));
}

void test_non_strict_find_block_by_global_name() {
  Document doc = parse("_a 1\ndata_foo\n", "<test>", false);
  const Block* g = doc.find_block(string_view("global_"));
  CHECK(g != 0);
  if (g) CHECK_EQ(sv_to_string(g->name()), std::string("global_"));
}

// ─── Explicit global_ header (CIF 1.1 reserved block type) ────────
// Separate from non-strict synthesis: `global_` appears as its own
// line and the tokenizer must recognize it as a block header, not a
// value. Case is insensitive. This is the form used by cctbx's
// monomer library (mon_lib_list.cif and friends).

void test_explicit_global_header_strict() {
  Document doc = parse("global_\n_a 1\n_b 2\n");
  CHECK_EQ(doc.size(), 1u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_a"))),
           std::string("1"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_b"))),
           std::string("2"));
}

void test_explicit_global_followed_by_data_block() {
  Document doc = parse("global_\n_a 1\ndata_foo\n_x 1\n");
  CHECK_EQ(doc.size(), 2u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_a"))),
           std::string("1"));
  CHECK_EQ(sv_to_string(doc[1].name()), std::string("foo"));
}

void test_explicit_global_case_insensitive() {
  Document doc = parse("GLOBAL_\n_a 1\n");
  CHECK_EQ(doc.size(), 1u);
  const Block* g = doc.find_block(string_view("global_"));
  CHECK(g != 0);
}

void test_monomer_library_header_with_explicit_global() {
  const char* src =
    "global_\n"
    "_lib_name         mon_lib\n"
    "_lib_version      4.11\n"
    "_lib_update       15/04/05\n"
    "# -----------------------------------------------\n"
    "data_comp_list\n"
    "loop_\n"
    "_chem_comp.id\n"
    "ALA\n";
  Document doc = parse(src, "<monomer>");  // strict=true — no relaxation needed
  CHECK_EQ(doc.size(), 2u);
  CHECK_EQ(sv_to_string(doc[0].name()), std::string("global_"));
  CHECK_EQ(sv_to_string(doc[0].find_value(string_view("_lib_version"))),
           std::string("4.11"));
  CHECK_EQ(sv_to_string(doc[1].name()), std::string("comp_list"));
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  test_strict_rejects_leading_pair();
  test_default_is_strict();
  test_strict_true_rejects_leading_pair();

  test_non_strict_leading_pair_then_data_block();
  test_non_strict_only_leading_pairs();
  test_non_strict_monomer_library_header_pattern();
  test_non_strict_leading_loop();
  test_non_strict_empty_input();
  test_non_strict_only_data_block_unchanged();
  test_non_strict_find_block_by_global_name();

  test_explicit_global_header_strict();
  test_explicit_global_followed_by_data_block();
  test_explicit_global_case_insensitive();
  test_monomer_library_header_with_explicit_global();
}

XCIF_TEST_MAIN()
