/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <stdio.h>
#include <cctbx/error.h>

namespace cctbx {

  error::error(std::string const& msg) throw()
  {
    msg_ = std::string("cctbx Error: ") + msg;
  }

  error::error(const char* file, long line, std::string const& msg,
               bool internal) throw()
  {
    const char *s = "";
    if (internal) s = " Internal";
    char buf[64];
    sprintf(buf, "%ld", line);
    msg_ =   std::string("cctbx") + s + " Error: "
              + file + "(" + buf + ")";
    if (msg.size()) msg_ += std::string(": ") + msg;
  }

  error::~error() throw() {}

  const char* error::what() const throw()
  {
     return msg_.c_str();
  }

  error_index::error_index(std::string const& msg) throw()
  : error(msg)
  {}

  error_index::~error_index() throw() {}

} // namespace cctbx
