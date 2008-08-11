#ifndef SCITBX_ARRAY_FAMILY_ERROR_H
#define SCITBX_ARRAY_FAMILY_ERROR_H

#include <stdexcept>

// FIXES for broken compilers
#include <boost/config.hpp>

namespace scitbx { namespace af {

  inline
  void throw_range_error()
  {
    throw std::range_error("scitbx array_family range error");
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ERROR_H
