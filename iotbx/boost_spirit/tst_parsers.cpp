#include <iotbx/boost_spirit/fortran_numerics.h>

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>

#include <iotbx/error.h>
#include <boost/optional.hpp>
#include <string>
#include <vector>
#include <iostream>

void exercise_fortran_int() {
  using namespace iotbx::boost_spirit;
  fortran_int_parser<int, /*Width=*/4> i4_p;
  std::vector<std::string> input;
  typedef boost::optional<int> opt_t;
  std::vector<opt_t> output;
  input.push_back("   1");  output.push_back(1);
  input.push_back("  1");  output.push_back(1);
  input.push_back(" 1");  output.push_back(1);
  input.push_back("1");  output.push_back(1);
  input.push_back("  1 ");  output.push_back(1);
  input.push_back("  12");  output.push_back(12);
  input.push_back(" 123");  output.push_back(123);
  input.push_back("1234");  output.push_back(1234);
  input.push_back(" 1 2");  output.push_back(12);
  input.push_back(" -1");  output.push_back(-1);
  input.push_back("-123");  output.push_back(-123);
  input.push_back("- 12");  output.push_back(-12);
  input.push_back("- 1 ");  output.push_back(-1);
  input.push_back("    ");  output.push_back(0);
  input.push_back("-1-2");  output.push_back(opt_t());
  input.push_back("1.024");  output.push_back(opt_t());
  for (int i=0; i < input.size(); ++i) {
    opt_t res;
    parse_info<std::string::iterator> info = parse(
      input[i].begin(), input[i].end(), i4_p[assign_a(res)]);
      /* Note: those parsers must be use only for character-level parsing
               since white spaces are significant */
    if (output[i]) {
      IOTBX_ASSERT(info.hit)(i)(input[i]);
      IOTBX_ASSERT(info.length == std::min(input[i].length(), std::size_t(4)))
                  (i)(input[i])(info.length);
      IOTBX_ASSERT(res == output[i])(i)(input[i])(*res);
    }
    else {
      IOTBX_ASSERT(!info.hit)(i)(input[i]);
    }
  }
  {
    std::vector<int> results;
    parse_info<> info = parse(
      " 1 2- 34", i4_p[push_back_a(results)] >> i4_p[push_back_a(results)]);
    IOTBX_ASSERT(results.size() == 2)(results.size());
    IOTBX_ASSERT(results[0] == 12)(results[0]);
    IOTBX_ASSERT(results[1] == -34)(results[1]);
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
  fortran_real_parser<double, /*Width=*/6, /*FracDigits=*/3> f63_p;
  std::vector<std::string> input;
  typedef boost::optional<double> opt_t;
  std::vector<opt_t> output;
  input.push_back("1."); output.push_back(1);
  input.push_back("12."); output.push_back(12);
  input.push_back("123."); output.push_back(123);
  input.push_back("12345."); output.push_back(12345);
  input.push_back("-1 3 ."); output.push_back(-13);
  input.push_back("1024"); output.push_back(1.024);
  input.push_back("102.4"); output.push_back(102.4);
  input.push_back(".1024"); output.push_back(0.1024);
  input.push_back("16+2"); output.push_back(1.6);
  input.push_back("16e+2"); output.push_back(1.6);
  input.push_back("16e2"); output.push_back(1.6);
  input.push_back("16.+2"); output.push_back(1600);
  input.push_back("16.e+2"); output.push_back(1600);
  input.push_back("16.e2"); output.push_back(1600);
  input.push_back("16-2"); output.push_back(16.e-5);
  input.push_back("16e-2"); output.push_back(16.e-5);
  input.push_back("16.-2"); output.push_back(0.16);
  input.push_back("16.e-2"); output.push_back(0.16);
  for (int i=0; i < input.size(); ++i) {
    opt_t res;
    parse_info<std::string::iterator> info = parse(
      input[i].begin(), input[i].end(),
      f63_p[assign_a(res)]);
    if (output[i]) {
      IOTBX_ASSERT(info.hit)(i)(input[i]);
      IOTBX_ASSERT(info.length == std::min(input[i].length(), std::size_t(6)))
                  (i)(input[i])(info.length);
      IOTBX_ASSERT(res == output[i])(i)(input[i])(*res);
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
      "0.1024-0.1289a",
      f63_p[push_back_a(results)] >> f63_p[push_back_a(results)]);
    IOTBX_ASSERT(info.hit);
    IOTBX_ASSERT(results[0] == 0.1024)(results[0]);
    IOTBX_ASSERT(results[1] == -0.128)(results[1]);
  }
  {
    std::vector<double> results;
    parse_info<> info = parse(
      "0.10240.2048a",
      f63_p[push_back_a(results)] >> f63_p[push_back_a(results)]);
    IOTBX_ASSERT(info.hit);
    IOTBX_ASSERT(results[0] == 0.1024)(results[0]);
    IOTBX_ASSERT(results[1] == 0.2048)(results[1]);
  }
  {
    std::vector<double> results;
    parse_info<> info = parse(
      "-16.-22",
      f63_p[push_back_a(results)] >> f63_p[push_back_a(results)]);
    IOTBX_ASSERT(info.hit);
    IOTBX_ASSERT(results[0] == -0.16)(results[0]);
    IOTBX_ASSERT(results[1] == 0.002)(results[1]);
  }
}

int main() {
  exercise_fortran_int();
  exercise_fortran_real_fixed();
  std::cout << "OK" << std::endl;
}
