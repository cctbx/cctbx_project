// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <stdio.h>
#include <cctbx/error.h>

namespace cctbx {

  error::error(const std::string& msg) throw()
  {
    message = std::string("cctbx Error: ") + msg;
  }

  const char* error::what() const throw() {
     return const_cast<char*>(message.c_str());
  }

  error::error(const char* file, long line, const std::string& msg) throw()
  {
    char buf[64];
    sprintf(buf, "%ld", line);
    message =   std::string("cctbx Internal Error: ")
              + file + "(" + buf + ")";
    if (msg.size()) message += std::string(": ") + msg;
  }

} // namespace cctbx
