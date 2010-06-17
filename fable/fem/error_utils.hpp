#ifndef FEM_ERROR_UTILS_HPP
#define FEM_ERROR_UTILS_HPP

#include <boost_adaptbx/error_utils.h>

#define FEM_THROW_UNHANDLED(message) \
  throw std::runtime_error( \
    boost_adaptbx::error_utils::file_and_line_as_string(__FILE__, __LINE__) \
    + ": UNHANDLED: " + message)

#endif // GUARD
