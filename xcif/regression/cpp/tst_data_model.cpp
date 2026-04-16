// cctbx_project/xcif/regression/cpp/tst_data_model.cpp
#include <xcif/data_model.h>
#include "test_utils.h"
#include <string>
#include <vector>

using namespace xcif;

// ─── §1  Document-level ────────────────────────────────────────────

void test_empty_document() {
  Document doc = parse("");
  CHECK_EQ(doc.size(), 0u);
}

void test_single_empty_block() {
  Document doc = parse("data_test");
  CHECK_EQ(doc.size(), 1u);
  CHECK_EQ(doc[0].name(), "test");
}

void test_multiple_blocks() {
  Document doc = parse("data_a\ndata_b\ndata_c\n");
  CHECK_EQ(doc.size(), 3u);
  CHECK_EQ(doc[0].name(), "a");
  CHECK_EQ(doc[1].name(), "b");
  CHECK_EQ(doc[2].name(), "c");
}

void test_find_block_by_name() {
  Document doc = parse("data_first\ndata_second\n");
  const Block* b = doc.find_block("second");
  CHECK(b != NULL);
  CHECK_EQ(b->name(), "second");
}

void test_find_block_case_insensitive() {
  Document doc = parse("data_MyBlock\n");
  CHECK(doc.find_block("myblock") != NULL);
  CHECK(doc.find_block("MYBLOCK") != NULL);
  CHECK(doc.find_block("MyBlock") != NULL);
}

void test_find_block_missing_returns_null() {
  Document doc = parse("data_a\n");
  CHECK(doc.find_block("nonexistent") == NULL);
}

// ─── §2  Tag-value pairs ───────────────────────────────────────────

void test_single_tag_value() {
  Document doc = parse("data_t\n_tag value\n");
  const Block& b = doc[0];
  CHECK(b.has_tag("_tag"));
  CHECK_EQ(b.find_value("_tag"), "value");
}

void test_multiple_tag_values() {
  Document doc = parse(
    "data_t\n"
    "_alpha 1\n"
    "_beta  2\n"
    "_gamma 5\n"
  );
  const Block& b = doc[0];
  CHECK_EQ(b.find_value("_alpha"), "1");
  CHECK_EQ(b.find_value("_beta"),  "2");
  CHECK_EQ(b.find_value("_gamma"), "5");
}

void test_tag_lookup_case_insensitive() {
  Document doc = parse("data_t\n_Some_Tag value\n");
  const Block& b = doc[0];
  CHECK_EQ(b.find_value("_some_tag"), "value");
  CHECK_EQ(b.find_value("_SOME_TAG"), "value");
}

void test_missing_tag_returns_empty() {
  Document doc = parse("data_t\n_x 1\n");
  CHECK_EQ(doc[0].find_value("_nope"), "");
}

void test_tag_value_separate_lines() {
  Document doc = parse(
    "data_t\n"
    "_tag\n"
    "value\n"
  );
  CHECK_EQ(doc[0].find_value("_tag"), "value");
}

void test_quoted_value() {
  Document doc = parse("data_t\n_q 'hello world'\n");
  CHECK_EQ(doc[0].find_value("_q"), "hello world");
}

void test_double_quoted_value() {
  Document doc = parse("data_t\n_q \"hello world\"\n");
  CHECK_EQ(doc[0].find_value("_q"), "hello world");
}

void test_semicolon_text_field() {
  Document doc = parse(
    "data_t\n_text\n;line one\nline two\n;\n"
  );
  CHECK_EQ(doc[0].find_value("_text"), "line one\nline two\n");
}

void test_dot_value() {
  Document doc = parse("data_t\n_x .\n");
  CHECK_EQ(doc[0].find_value("_x"), ".");
}

void test_question_mark_value() {
  Document doc = parse("data_t\n_x ?\n");
  CHECK_EQ(doc[0].find_value("_x"), "?");
}

// ─── §3  Loops ─────────────────────────────────────────────────────

void test_simple_loop() {
  Document doc = parse(
    "data_t\n"
    "loop_\n"
    "_a\n_b\n"
    "1 2\n"
    "3 4\n"
  );
  const Loop* lp = doc[0].find_loop("_a");
  CHECK(lp != NULL);
  CHECK_EQ(lp->width(), 2u);
  CHECK_EQ(lp->length(), 2u);
}

void test_loop_tags() {
  Document doc = parse(
    "data_t\n"
    "loop_\n_x\n_y\n_z\n"
    "1 2 3\n"
  );
  const Loop* lp = doc[0].find_loop("_x");
  CHECK(lp != NULL);
  const std::vector<string_view>& tags = lp->tags();
  CHECK_EQ(tags.size(), 3u);
  CHECK_EQ(tags[0], "_x");
  CHECK_EQ(tags[1], "_y");
  CHECK_EQ(tags[2], "_z");
}

void test_loop_value_access() {
  Document doc = parse(
    "data_t\n"
    "loop_\n_a\n_b\n"
    "x1 y1\n"
    "x2 y2\n"
  );
  const Loop* lp = doc[0].find_loop("_a");
  CHECK(lp != NULL);
  CHECK_EQ(lp->value(0, 0), "x1");
  CHECK_EQ(lp->value(0, 1), "y1");
  CHECK_EQ(lp->value(1, 0), "x2");
  CHECK_EQ(lp->value(1, 1), "y2");
}

void test_loop_column_by_tag() {
  Document doc = parse(
    "data_t\n"
    "loop_\n_id\n_val\n"
    "1 a\n2 b\n3 c\n"
  );
  const Loop* lp = doc[0].find_loop("_id");
  CHECK(lp != NULL);
  std::vector<string_view> col = lp->column("_val");
  CHECK_EQ(col.size(), 3u);
  CHECK_EQ(col[0], "a");
  CHECK_EQ(col[1], "b");
  CHECK_EQ(col[2], "c");
}

void test_loop_column_case_insensitive() {
  Document doc = parse(
    "data_t\nloop_\n_Tag\n1\n2\n"
  );
  const Loop* lp = doc[0].find_loop("_tag");
  CHECK(lp != NULL);
  CHECK_EQ(lp->column("_TAG").size(), 2u);
}

void test_loop_has_tag() {
  Document doc = parse(
    "data_t\nloop_\n_a\n_b\n1 2\n"
  );
  const Loop* lp = doc[0].find_loop("_a");
  CHECK(lp != NULL);
  CHECK(lp->has_tag("_a"));
  CHECK(lp->has_tag("_b"));
  CHECK(!lp->has_tag("_c"));
}

void test_find_loop_returns_null_for_scalar() {
  Document doc = parse("data_t\n_x 1\n");
  CHECK(doc[0].find_loop("_x") == NULL);
}

void test_find_loop_returns_null_for_missing() {
  Document doc = parse("data_t\n_x 1\n");
  CHECK(doc[0].find_loop("_nope") == NULL);
}

void test_multiple_loops() {
  Document doc = parse(
    "data_t\n"
    "loop_\n_a\n1\n2\n"
    "loop_\n_b\n3\n4\n"
  );
  const Block& b = doc[0];
  CHECK(b.find_loop("_a") != NULL);
  CHECK(b.find_loop("_b") != NULL);
  CHECK(b.find_loop("_a") != b.find_loop("_b"));
  CHECK_EQ(b.loops().size(), 2u);
}

void test_single_column_loop() {
  Document doc = parse(
    "data_t\nloop_\n_only\na\nb\nc\n"
  );
  const Loop* lp = doc[0].find_loop("_only");
  CHECK(lp != NULL);
  CHECK_EQ(lp->width(), 1u);
  CHECK_EQ(lp->length(), 3u);
}

void test_loop_terminated_by_new_block() {
  Document doc = parse(
    "data_a\nloop_\n_x\n1\n2\n"
    "data_b\n_y 3\n"
  );
  CHECK_EQ(doc.size(), 2u);
  CHECK_EQ(doc[0].find_loop("_x")->length(), 2u);
  CHECK_EQ(doc[1].find_value("_y"), "3");
}

// ─── §4  Mixed content ─────────────────────────────────────────────

void test_tags_and_loops_together() {
  Document doc = parse(
    "data_t\n"
    "_scalar 42\n"
    "loop_\n_col\na\nb\n"
    "_another val\n"
  );
  const Block& b = doc[0];
  CHECK_EQ(b.find_value("_scalar"), "42");
  CHECK_EQ(b.find_value("_another"), "val");
  CHECK(b.find_loop("_col") != NULL);
  CHECK_EQ(b.find_loop("_col")->length(), 2u);
}

void test_loop_with_quoted_values() {
  Document doc = parse(
    "data_t\nloop_\n_name\n'Alice Bob'\n\"Carol\"\n"
  );
  const Loop* lp = doc[0].find_loop("_name");
  CHECK(lp != NULL);
  CHECK_EQ(lp->value(0, 0), "Alice Bob");
  CHECK_EQ(lp->value(1, 0), "Carol");
}

// ─── §5  Save frames ───────────────────────────────────────────────

void test_save_frame() {
  Document doc = parse(
    "data_dict\n"
    "save_my_frame\n"
    "_tag val\n"
    "save_\n"
  );
  const Block& b = doc[0];
  const Block* sf = b.find_save_frame("my_frame");
  CHECK(sf != NULL);
  CHECK_EQ(sf->name(), "my_frame");
  CHECK_EQ(sf->find_value("_tag"), "val");
}

void test_save_frame_case_insensitive() {
  Document doc = parse(
    "data_d\nsave_Frame\n_x 1\nsave_\n"
  );
  CHECK(doc[0].find_save_frame("frame") != NULL);
  CHECK(doc[0].find_save_frame("FRAME") != NULL);
}

void test_save_frame_with_loop() {
  Document doc = parse(
    "data_d\n"
    "save_s\n"
    "loop_\n_a\n1\n2\n"
    "save_\n"
  );
  const Block* sf = doc[0].find_save_frame("s");
  CHECK(sf != NULL);
  const Loop* lp = sf->find_loop("_a");
  CHECK(lp != NULL);
  CHECK_EQ(lp->length(), 2u);
}

void test_multiple_save_frames() {
  Document doc = parse(
    "data_d\n"
    "save_one\n_x 1\nsave_\n"
    "save_two\n_y 2\nsave_\n"
  );
  CHECK(doc[0].find_save_frame("one") != NULL);
  CHECK(doc[0].find_save_frame("two") != NULL);
}

// ─── §6  Misc ───────────────────────────────────────────────────────

void test_block_name_preserves_original_case() {
  Document doc = parse("data_MyProtein\n");
  CHECK_EQ(doc[0].name(), "MyProtein");
}

void test_block_source_order() {
  Document doc = parse(
    "data_t\n"
    "_a 1\n"
    "loop_ _x _y 1 2\n"
    "_b 2\n",
    "<test>");
  const Block& blk = doc[0];
  const auto& order = blk.source_order();
  CHECK_EQ(order.size(), 3u);
  CHECK_EQ(order[0].first, Block::ENTRY_PAIR);
  CHECK_EQ(order[0].second, 0u);
  CHECK_EQ(order[1].first, Block::ENTRY_LOOP);
  CHECK_EQ(order[1].second, 0u);
  CHECK_EQ(order[2].first, Block::ENTRY_PAIR);
  CHECK_EQ(order[2].second, 1u);
}

void test_block_source_order_with_save_frame() {
  Document doc = parse(
    "data_t\n"
    "_a 1\n"
    "save_s1\n"
    "_inside 2\n"
    "save_\n"
    "loop_ _x _y 1 2\n",
    "<test>");
  const Block& blk = doc[0];
  const auto& order = blk.source_order();
  CHECK_EQ(order.size(), 3u);
  CHECK_EQ(order[0].first, Block::ENTRY_PAIR);
  CHECK_EQ(order[0].second, 0u);
  CHECK_EQ(order[1].first, Block::ENTRY_SAVE_FRAME);
  CHECK_EQ(order[1].second, 0u);
  CHECK_EQ(order[2].first, Block::ENTRY_LOOP);
  CHECK_EQ(order[2].second, 0u);
  // Save frame's own source_order should contain its pair item.
  const Block& sf = blk.save_frames()[0];
  const auto& sf_order = sf.source_order();
  CHECK_EQ(sf_order.size(), 1u);
  CHECK_EQ(sf_order[0].first, Block::ENTRY_PAIR);
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  // §1 Document-level
  test_empty_document();
  test_single_empty_block();
  test_multiple_blocks();
  test_find_block_by_name();
  test_find_block_case_insensitive();
  test_find_block_missing_returns_null();

  // §2 Tag-value pairs
  test_single_tag_value();
  test_multiple_tag_values();
  test_tag_lookup_case_insensitive();
  test_missing_tag_returns_empty();
  test_tag_value_separate_lines();
  test_quoted_value();
  test_double_quoted_value();
  test_semicolon_text_field();
  test_dot_value();
  test_question_mark_value();

  // §3 Loops
  test_simple_loop();
  test_loop_tags();
  test_loop_value_access();
  test_loop_column_by_tag();
  test_loop_column_case_insensitive();
  test_loop_has_tag();
  test_find_loop_returns_null_for_scalar();
  test_find_loop_returns_null_for_missing();
  test_multiple_loops();
  test_single_column_loop();
  test_loop_terminated_by_new_block();

  // §4 Mixed content
  test_tags_and_loops_together();
  test_loop_with_quoted_values();

  // §5 Save frames
  test_save_frame();
  test_save_frame_case_insensitive();
  test_save_frame_with_loop();
  test_multiple_save_frames();

  // §6 Misc
  test_block_name_preserves_original_case();
  test_block_source_order();
  test_block_source_order_with_save_frame();
}

XCIF_TEST_MAIN()
