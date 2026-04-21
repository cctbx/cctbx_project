// cctbx_project/xcif/regression/cpp/tst_numeric.cpp
#include <xcif/numeric.h>
#include "test_utils.h"
#include <cmath>
#include <limits>

using xcif::string_view;
using xcif::as_double;
using xcif::as_int;
using xcif::as_double_with_su;
using xcif::is_unknown;
using xcif::is_inapplicable;
using xcif::is_null;

static bool approx_eq(double a, double b, double eps = 1e-12) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::fabs(a - b) < eps;
}

// ─── §1  as_double — common path (no SU) ───────────────────────────

void test_double_integer() {
  CHECK(approx_eq(as_double("42"), 42.0));
}

void test_double_simple_float() {
  CHECK(approx_eq(as_double("3.14"), 3.14));
}

void test_double_negative() {
  CHECK(approx_eq(as_double("-1.5"), -1.5));
}

void test_double_leading_plus() {
  CHECK(approx_eq(as_double("+1.5"), 1.5));
}

void test_double_zero() {
  CHECK(approx_eq(as_double("0"), 0.0));
}

void test_double_zero_float() {
  CHECK(approx_eq(as_double("0.0"), 0.0));
}

void test_double_scientific_lower() {
  CHECK(approx_eq(as_double("1.5e3"), 1500.0));
}

void test_double_scientific_upper() {
  CHECK(approx_eq(as_double("1.5E3"), 1500.0));
}

void test_double_scientific_neg_exp() {
  CHECK(approx_eq(as_double("2.5e-4"), 2.5e-4));
}

void test_double_scientific_plus_exp() {
  CHECK(approx_eq(as_double("1.5e+3"), 1500.0));
}

void test_double_large_value() {
  CHECK(approx_eq(as_double("1.23456789e+20"), 1.23456789e+20, 1e8));
}

void test_double_small_value() {
  CHECK(approx_eq(as_double("1.23e-30"), 1.23e-30, 1e-42));
}

void test_double_no_leading_digit() {
  CHECK(approx_eq(as_double(".5"), 0.5));
}

void test_double_negative_no_leading_digit() {
  CHECK(approx_eq(as_double("-.5"), -0.5));
}

// ─── §2  as_double — SU stripping (rare path) ──────────────────────

void test_double_strips_simple_su() {
  CHECK(approx_eq(as_double("1.542(3)"), 1.542));
}

void test_double_strips_multidigit_su() {
  CHECK(approx_eq(as_double("23.45(12)"), 23.45));
}

void test_double_strips_integer_su() {
  CHECK(approx_eq(as_double("100(5)"), 100.0));
}

void test_double_strips_su_many_decimals() {
  CHECK(approx_eq(as_double("1.54321(5)"), 1.54321));
}

// ─── §3  as_double — special values ─────────────────────────────────

void test_double_unknown_is_nan() {
  CHECK(std::isnan(as_double(".")));
}

void test_double_inapplicable_is_nan() {
  CHECK(std::isnan(as_double("?")));
}

void test_double_non_numeric_is_nan() {
  CHECK(std::isnan(as_double("hello")));
}

// ─── §4  as_int ─────────────────────────────────────────────────────

void test_int_simple() {
  CHECK_EQ(as_int("42"), 42);
}

void test_int_negative() {
  CHECK_EQ(as_int("-7"), -7);
}

void test_int_zero() {
  CHECK_EQ(as_int("0"), 0);
}

void test_int_large() {
  CHECK_EQ(as_int("123456"), 123456);
}

void test_int_leading_plus() {
  CHECK_EQ(as_int("+5"), 5);
}

void test_int_max() {
  CHECK_EQ(as_int("2147483647"), std::numeric_limits<int>::max());
}

void test_int_min() {
  CHECK_EQ(as_int("-2147483648"), std::numeric_limits<int>::min());
}

void test_int_overflow_pos() {
  bool caught = false;
  try { as_int("2147483648"); }
  catch (const std::overflow_error&) { caught = true; }
  CHECK(caught);
}

void test_int_overflow_neg() {
  bool caught = false;
  try { as_int("-2147483649"); }
  catch (const std::overflow_error&) { caught = true; }
  CHECK(caught);
}

// ─── §5  as_double_with_su ──────────────────────────────────────────

void test_su_none() {
  std::pair<double, double> r = as_double_with_su("1.542");
  CHECK(approx_eq(r.first, 1.542));
  CHECK(approx_eq(r.second, 0.0));
}

void test_su_single_digit() {
  std::pair<double, double> r = as_double_with_su("1.542(3)");
  CHECK(approx_eq(r.first, 1.542));
  CHECK(approx_eq(r.second, 0.003));
}

void test_su_multidigit() {
  std::pair<double, double> r = as_double_with_su("23.45(12)");
  CHECK(approx_eq(r.first, 23.45));
  CHECK(approx_eq(r.second, 0.12));
}

void test_su_integer() {
  std::pair<double, double> r = as_double_with_su("100(5)");
  CHECK(approx_eq(r.first, 100.0));
  CHECK(approx_eq(r.second, 5.0));
}

void test_su_many_decimals() {
  std::pair<double, double> r = as_double_with_su("1.54321(5)");
  CHECK(approx_eq(r.first, 1.54321));
  CHECK(approx_eq(r.second, 0.00005));
}

void test_su_single_decimal() {
  std::pair<double, double> r = as_double_with_su("1.0(1)");
  CHECK(approx_eq(r.first, 1.0));
  CHECK(approx_eq(r.second, 0.1));
}

void test_su_with_exponent() {
  // "1.5e+3(2)": 1 decimal digit before 'e', so su = 2/10 = 0.2; value = 1500
  std::pair<double, double> r = as_double_with_su("1.5e+3(2)");
  CHECK(approx_eq(r.first, 1500.0));
  CHECK(approx_eq(r.second, 0.2));
}

void test_su_unknown_returns_nan() {
  std::pair<double, double> r = as_double_with_su(".");
  CHECK(std::isnan(r.first));
  CHECK(approx_eq(r.second, 0.0));
}

// ─── §5b as_double / as_double_with_su — boundary / malformed input ──

void test_double_long_value_nan() {
  // Values with ≥32 chars before '(' exceed the stack buffer → NaN
  CHECK(std::isnan(as_double("1.23456789012345678901234567890123")));
}

void test_su_missing_close_paren() {
  // "1.5(" has no closing ')' → su = 0.0 (graceful)
  std::pair<double, double> r = as_double_with_su("1.5(");
  CHECK(approx_eq(r.first, 1.5));
  CHECK(approx_eq(r.second, 0.0));
}

void test_su_nonnumeric_in_su_content() {
  // Non-digit chars in su content are skipped; digits are accumulated.
  // "1.5(3a5)": su_int = 35, decimals = 1, su = 3.5
  std::pair<double, double> r = as_double_with_su("1.5(3a5)");
  CHECK(approx_eq(r.first, 1.5));
  CHECK(approx_eq(r.second, 3.5));
}

// ─── §6  Predicates ─────────────────────────────────────────────────

void test_is_unknown_dot() {
  CHECK(is_unknown("."));
}

void test_is_unknown_other() {
  CHECK(!is_unknown("x"));
  CHECK(!is_unknown("?"));
  CHECK(!is_unknown(".."));
}

void test_is_inapplicable_question() {
  CHECK(is_inapplicable("?"));
}

void test_is_inapplicable_other() {
  CHECK(!is_inapplicable("x"));
  CHECK(!is_inapplicable("."));
}

void test_is_null_both() {
  CHECK(is_null("."));
  CHECK(is_null("?"));
  CHECK(!is_null("x"));
  CHECK(!is_null("1.5"));
}

// ─── run_all_tests ─────────────────────────────────────────────────

void run_all_tests() {
  // §1 as_double — common path
  test_double_integer();
  test_double_simple_float();
  test_double_negative();
  test_double_leading_plus();
  test_double_zero();
  test_double_zero_float();
  test_double_scientific_lower();
  test_double_scientific_upper();
  test_double_scientific_neg_exp();
  test_double_scientific_plus_exp();
  test_double_large_value();
  test_double_small_value();
  test_double_no_leading_digit();
  test_double_negative_no_leading_digit();

  // §2 as_double — SU stripping
  test_double_strips_simple_su();
  test_double_strips_multidigit_su();
  test_double_strips_integer_su();
  test_double_strips_su_many_decimals();

  // §3 as_double — special values
  test_double_unknown_is_nan();
  test_double_inapplicable_is_nan();
  test_double_non_numeric_is_nan();

  // §4 as_int
  test_int_simple();
  test_int_negative();
  test_int_zero();
  test_int_large();
  test_int_leading_plus();
  test_int_max();
  test_int_min();
  test_int_overflow_pos();
  test_int_overflow_neg();

  // §5 as_double_with_su
  test_su_none();
  test_su_single_digit();
  test_su_multidigit();
  test_su_integer();
  test_su_many_decimals();
  test_su_single_decimal();
  test_su_with_exponent();
  test_su_unknown_returns_nan();

  // §5b as_double/as_double_with_su — boundary / malformed
  test_double_long_value_nan();
  test_su_missing_close_paren();
  test_su_nonnumeric_in_su_content();

  // §6 Predicates
  test_is_unknown_dot();
  test_is_unknown_other();
  test_is_inapplicable_question();
  test_is_inapplicable_other();
  test_is_null_both();
}

XCIF_TEST_MAIN()
