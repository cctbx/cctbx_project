// $Id$

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
