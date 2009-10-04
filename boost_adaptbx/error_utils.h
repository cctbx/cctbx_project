#ifndef BOOST_ADAPTBX_ERROR_UTILS_H
#define BOOST_ADAPTBX_ERROR_UTILS_H

#include <sstream>
#include <string>
#include <stdexcept>

namespace boost_adaptbx { namespace error_utils {

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

}} // boost_adaptbx::error_utils

#define ASSERTBX(condition) \
  if (!(condition)) { \
    throw std::runtime_error( \
      boost_adaptbx::error_utils::file_and_line_as_string(__FILE__, __LINE__) \
      + ": ASSERT(" #condition ") failure."); \
  }

#define BOOST_ADAPTBX_UNREACHABLE_ERROR() \
  std::runtime_error( \
    "Control flow passes through branch that should be unreachable: " \
    + boost_adaptbx::error_utils::file_and_line_as_string(__FILE__, __LINE__))

#endif // GUARD
