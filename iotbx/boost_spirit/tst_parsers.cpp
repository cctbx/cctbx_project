#include <iotbx/boost_spirit/fortran_int.h>
#include <iotbx/boost_spirit/fortran_real_fixed.h>
#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <iotbx/error.h>
#include <boost/optional.hpp>
#include <string>
#include <vector>
#include <iostream>

void exercise_fortran_int() {
  using namespace iotbx::boost_spirit;
  fortran_int_parser<4> i4_p;
  std::vector<std::string> input;
  typedef boost::optional<int> opt_t;
  std::vector<opt_t> output;
  input.push_back("   1");  output.push_back(1);
  input.push_back("  12");  output.push_back(12);
  input.push_back(" 123");  output.push_back(123);
  input.push_back("1234");  output.push_back(1234);
  input.push_back("    1");  output.push_back(opt_t());
  input.push_back("   12");  output.push_back(1);
  input.push_back("  123");  output.push_back(12);
  input.push_back(" 1234");  output.push_back(123);
  input.push_back("12345"); output.push_back(1234);
  input.push_back("1234a"); output.push_back(1234);
  input.push_back("  -1");  output.push_back(-1);
  input.push_back(" -12");  output.push_back(-12);
  input.push_back("-123");  output.push_back(-123);
  input.push_back("-1234"); output.push_back(-123);
  input.push_back("-123a"); output.push_back(-123);
  input.push_back("    ");  output.push_back(opt_t());
  input.push_back("  a ");  output.push_back(opt_t());
  input.push_back("--12");  output.push_back(opt_t());
  input.push_back(" 1 2");  output.push_back(opt_t());
  input.push_back("   -");  output.push_back(opt_t());
  input.push_back("123");  output.push_back(opt_t());
  input.push_back("123");  output.push_back(opt_t());
  for (int i=0; i < input.size(); ++i) {
    opt_t res;
    parse_info<std::string::iterator> info = parse(
      input[i].begin(), input[i].end(), i4_p[assign_a(res)]);
      /* Note: those parsers must be use only for character-level parsing
               since white spaces are significant */
    if (output[i]) {
      IOTBX_ASSERT(info.hit)(i)(input[i]);
      IOTBX_ASSERT(info.length == 4)(i)(input[i])(info.length);
      IOTBX_ASSERT(res == output[i])(i)(input[i])(res);
    }
    else {
      IOTBX_ASSERT(!info.hit)(i)(input[i]);
    }
  }
  {
    std::vector<int> results;
    parse_info<> info = parse(
      "  12  34", i4_p[push_back_a(results)] >> i4_p[push_back_a(results)]);
    IOTBX_ASSERT(results.size() == 2)(results.size());
    IOTBX_ASSERT(results[0] == 12)(results[0]);
    IOTBX_ASSERT(results[1] == 34)(results[1]);
  }
  {
    std::vector<int> results;
    parse_info<> info = parse(
      "1234-567", i4_p[push_back_a(results)] >> i4_p[push_back_a(results)]);
    IOTBX_ASSERT(results.size() == 2)(results.size());
    IOTBX_ASSERT(results[0] == 1234)(results[0]);
    IOTBX_ASSERT(results[1] == -567)(results[1]);
  }
}

void exercise_fortran_real_fixed() {
  using namespace iotbx::boost_spirit;
  fortran_real_fixed</*Width=*/6, /*FracDigits=*/3> f63_p;
  std::vector<std::string> input;
  typedef boost::optional<int> opt_t;
  std::vector<opt_t> output;
  input.push_back(" 1.024");  output.push_back(1.024);
  input.push_back("21.024");  output.push_back(21.024);
  input.push_back("-1.024");  output.push_back(-1.024);
  input.push_back(" 1.0245");  output.push_back(1.024);
  input.push_back(" 1.024a");  output.push_back(1.024);
  input.push_back("  1.024");  output.push_back(opt_t());
  input.push_back("123.024");  output.push_back(opt_t());
  input.push_back("123.24");  output.push_back(opt_t());
  for (int i=0; i < input.size(); ++i) {
    opt_t res;
    parse_info<std::string::iterator> info = parse(
      input[i].begin(), input[i].end(),
      f63_p[assign_a(res)]);
    if (output[i]) {
      IOTBX_ASSERT(info.hit)(i)(input[i]);
      IOTBX_ASSERT(info.length == 6)(i)(input[i])(info.length);
      IOTBX_ASSERT(res == output[i])(i)(input[i])(res);
    }
    else {
      IOTBX_ASSERT(!info.hit)(i)(input[i]);
    }
  }
  {
    std::vector<double> results;
    parse_info<> info = parse(
      " 1.024 2.048",
      f63_p[push_back_a(results)] >> f63_p[push_back_a(results)]);
    IOTBX_ASSERT(results.size() == 2)(results.size());
    IOTBX_ASSERT(results[0] == 1.024)(results[0]);
    IOTBX_ASSERT(results[1] == 2.048)(results[1]);
  }
  {
    std::vector<double> results;
    parse_info<> info = parse(
      "11.024-2.048a",
      f63_p[push_back_a(results)] >> f63_p[push_back_a(results)]);
    IOTBX_ASSERT(results.size() == 2)(results.size());
    IOTBX_ASSERT(results[0] == 11.024)(results[0]);
    IOTBX_ASSERT(results[1] == -2.048)(results[1]);
  }
  {
    std::vector<double> results;
    parse_info<> info = parse(
      "11.02422.048a",
      f63_p[push_back_a(results)] >> f63_p[push_back_a(results)]);
    IOTBX_ASSERT(results.size() == 2)(results.size());
    IOTBX_ASSERT(results[0] == 11.024)(results[0]);
    IOTBX_ASSERT(results[1] == 22.048)(results[1]);
  }
  {
    parse_info<> info = parse("11.024222.048 ", f63_p >> f63_p);
    IOTBX_ASSERT(!info.hit);
  }
}

int main() {
  exercise_fortran_int();
  exercise_fortran_real_fixed();
  std::cout << "OK" << std::endl;
}
