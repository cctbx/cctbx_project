// cctbx_project/xcif/regression/cpp/tst_string_view.cpp
#include <xcif/string_view.h>
#include "test_utils.h"
#include <string>
#include <sstream>

using xcif::string_view;

// ─── §1  Construction ──────────────────────────────────────────────

void test_default_construction() {
  string_view sv;
  CHECK(sv.empty());
  CHECK_EQ(sv.size(), 0u);
}

void test_from_ptr_and_length() {
  const char* buf = "hello world";
  string_view sv(buf, 5);
  CHECK_EQ(sv.size(), 5u);
  CHECK_EQ(sv, "hello");
}

void test_from_c_string() {
  string_view sv("test");
  CHECK_EQ(sv.size(), 4u);
  CHECK_EQ(sv, "test");
}

void test_from_null_c_string() {
  string_view sv(static_cast<const char*>(0));
  CHECK(sv.empty());
  CHECK_EQ(sv.size(), 0u);
}

void test_from_std_string() {
  std::string s("hello");
  string_view sv(s);
  CHECK_EQ(sv.size(), 5u);
  CHECK_EQ(sv, "hello");
}

// ─── §2  Element access ────────────────────────────────────────────

void test_subscript() {
  string_view sv("abc");
  CHECK_EQ(sv[0], 'a');
  CHECK_EQ(sv[1], 'b');
  CHECK_EQ(sv[2], 'c');
}

void test_front_and_back() {
  string_view sv("xyz");
  CHECK_EQ(sv.front(), 'x');
  CHECK_EQ(sv.back(), 'z');
}

void test_data_points_to_original() {
  const char* buf = "original";
  string_view sv(buf, 8);
  CHECK_EQ(sv.data(), buf);
}

// ─── §3  Comparison ────────────────────────────────────────────────

void test_equal_views() {
  string_view a("abc");
  string_view b("abc");
  CHECK(a == b);
  CHECK(!(a != b));
}

void test_unequal_views() {
  string_view a("abc");
  string_view b("xyz");
  CHECK(a != b);
  CHECK(!(a == b));
}

void test_unequal_lengths() {
  string_view a("abc");
  string_view b("ab");
  CHECK(a != b);
}

void test_compare_with_c_string() {
  string_view sv("hello");
  CHECK(sv == "hello");
  CHECK("hello" == sv);
  CHECK(sv != "world");
  CHECK("world" != sv);
}

void test_compare_with_std_string() {
  string_view sv("hello");
  std::string s("hello");
  CHECK(sv == s);
  CHECK(s == sv);
  CHECK(sv != std::string("other"));
}

void test_less_than() {
  CHECK(string_view("abc") < string_view("abd"));
  CHECK(string_view("ab") < string_view("abc"));
  CHECK(!(string_view("abc") < string_view("abc")));
}

// ─── §4  Operations ────────────────────────────────────────────────

void test_substr() {
  string_view sv("hello world");
  string_view sub = sv.substr(6, 5);
  CHECK_EQ(sub, "world");
}

void test_substr_to_end() {
  string_view sv("hello world");
  string_view sub = sv.substr(6);
  CHECK_EQ(sub, "world");
}

void test_substr_clamps() {
  string_view sv("abc");
  string_view sub = sv.substr(1, 100);
  CHECK_EQ(sub, "bc");
}

void test_find_char() {
  string_view sv("a.b.c");
  CHECK_EQ(sv.find('.'), 1u);
  CHECK_EQ(sv.find('.', 2), 3u);
  CHECK_EQ(sv.find('z'), string_view::npos);
}

void test_remove_prefix() {
  string_view sv("hello");
  sv.remove_prefix(2);
  CHECK_EQ(sv, "llo");
}

void test_remove_suffix() {
  string_view sv("hello");
  sv.remove_suffix(2);
  CHECK_EQ(sv, "hel");
}

// ─── §5  Conversion ────────────────────────────────────────────────

void test_explicit_string_conversion() {
  string_view sv("test");
  std::string s(sv);  // explicit via direct initialization
  CHECK_EQ(s, "test");
}

void test_view_into_buffer_does_not_copy() {
  const char buf[] = "buffer content";
  string_view sv(buf, 6);
  CHECK_EQ(sv.data(), buf);
  CHECK_EQ(sv, "buffer");
}

// ─── §6  Stream output ─────────────────────────────────────────────

void test_stream_output() {
  string_view sv("hello");
  std::ostringstream oss;
  oss << sv;
  CHECK_EQ(oss.str(), "hello");
}

void test_stream_output_empty() {
  string_view sv;
  std::ostringstream oss;
  oss << sv;
  CHECK_EQ(oss.str(), "");
}

// ─── §7  Iterator ───────────────────────────────────────────────────

void test_begin_end() {
  string_view sv("abc");
  std::string built;
  for (string_view::const_iterator it = sv.begin(); it != sv.end(); ++it) {
    built += *it;
  }
  CHECK_EQ(built, "abc");
}

void test_range_for() {
  string_view sv("xyz");
  std::string built;
  for (char c : sv) {
    built += c;
  }
  CHECK_EQ(built, "xyz");
}

// ─── §8  Empty and zero-length views ────────────────────────────────

void test_empty_string_view() {
  string_view sv("", 0);
  CHECK(sv.empty());
  CHECK_EQ(sv, "");
}

void test_empty_equals_empty() {
  string_view a;
  string_view b("", 0);
  CHECK(a == b);
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  // §1 Construction
  test_default_construction();
  test_from_ptr_and_length();
  test_from_c_string();
  test_from_null_c_string();
  test_from_std_string();

  // §2 Element access
  test_subscript();
  test_front_and_back();
  test_data_points_to_original();

  // §3 Comparison
  test_equal_views();
  test_unequal_views();
  test_unequal_lengths();
  test_compare_with_c_string();
  test_compare_with_std_string();
  test_less_than();

  // §4 Operations
  test_substr();
  test_substr_to_end();
  test_substr_clamps();
  test_find_char();
  test_remove_prefix();
  test_remove_suffix();

  // §5 Conversion
  test_explicit_string_conversion();
  test_view_into_buffer_does_not_copy();

  // §6 Stream output
  test_stream_output();
  test_stream_output_empty();

  // §7 Iterator
  test_begin_end();
  test_range_for();

  // §8 Edge cases
  test_empty_string_view();
  test_empty_equals_empty();
}

XCIF_TEST_MAIN()
