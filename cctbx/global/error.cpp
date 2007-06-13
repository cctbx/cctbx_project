#include <cctbx/error.h>

#undef CCTBX_ASSERT_HACK_A
#undef CCTBX_ASSERT_HACK_B

namespace cctbx {

  error::error(std::string const& msg) throw()
    : CCTBX_ASSERT_HACK_A(*this), CCTBX_ASSERT_HACK_B(*this)
  {
    msg_ = std::string("cctbx Error: ") + msg;
  }

  error::error(const char* file, long line, std::string const& msg,
               bool internal) throw()
    : CCTBX_ASSERT_HACK_A(*this), CCTBX_ASSERT_HACK_B(*this)
  {
     std::ostringstream o;
     o << "cctbx";
     if (internal) o << " Internal";
     o << "Error: " << file << "(" << line << ")";
     if (msg.size()) o << ": " << msg;
     msg_ = o.str();
//     const char *s = "";
//    if (internal) s = " Internal";
//    char buf[64];
//    sprintf(buf, "%ld", line);
//    msg_ =   std::string("cctbx") + s + " Error: "
//              + file + "(" + buf + ")";
//    if (msg.size()) msg_ += std::string(": ") + msg;
  }

  error::~error() throw() {}

  const char* error::what() const throw()
  {
    std::string result = msg_;
    if (values_.length() > 0) result += "\n" + values_;
    return result.c_str();
  }

  error_index::error_index(std::string const& msg) throw()
  : error(msg)
  {}

  error_index::~error_index() throw() {}

} // namespace cctbx
