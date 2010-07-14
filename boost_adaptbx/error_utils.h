#ifndef BOOST_ADAPTBX_ERROR_UTILS_H
#define BOOST_ADAPTBX_ERROR_UTILS_H

#include <boost_adaptbx/libc_backtrace.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

#define BOOSTBX_CHECK_POINT \
  std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl
#define BOOSTBX_CHECK_POINT_MSG(msg) \
  std::cout << msg << " @ " __FILE__ << "(" << __LINE__ << ")" << std::endl
#define BOOSTBX_EXAMINE(A) \
  std::cout << "variable " << (#A) << ": " << A << std::endl

namespace boost_adaptbx { namespace error_utils {

  //! Debugging aid: add call to code, recompile, run in gdb, use "where"
  inline
  int
  segfault_if(
    bool condition)
  {
    if (condition) {
      int* ptr = 0;
      return *ptr;
    }
    return 0;
  }

  inline
  std::string
  file_and_line_as_string(
    const char* file,
    long line)
  {
    std::ostringstream o;
    o << file << "(" << line << ")";
    return o.str();
  }

  inline
  std::string
  file_and_line_as_string_with_backtrace(
    const char* file,
    long line)
  {
    std::ostringstream o;
    boost_adaptbx::libc_backtrace::show_if_possible(o, 1);
    o << file << "(" << line << ")";
    return o.str();
  }

}} // boost_adaptbx::error_utils

#define ASSERTBX(condition) \
  if (!(condition)) { \
    throw std::runtime_error( \
      boost_adaptbx::error_utils::file_and_line_as_string( \
        __FILE__, __LINE__) \
      + ": ASSERT(" #condition ") failure."); \
  }

#define BOOST_ADAPTBX_UNREACHABLE_ERROR() \
  std::runtime_error( \
    "Control flow passes through branch that should be unreachable: " \
    + boost_adaptbx::error_utils::file_and_line_as_string(__FILE__, __LINE__))

#define BOOST_ADAPTBX_NOT_IMPLEMENTED() \
  std::runtime_error( \
    "Not implemented: " \
    + boost_adaptbx::error_utils::file_and_line_as_string(__FILE__, __LINE__))

#endif // GUARD
