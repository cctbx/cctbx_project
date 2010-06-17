#ifndef FEM_STOP_HPP
#define FEM_STOP_HPP

#include <fem/size_t.hpp>
#include <boost/shared_array.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fem {

  struct stop_info : public std::exception
  {
    std::string source_file;
    size_t source_line;
    ssize_t digits;
    std::string message;

    protected:
      mutable boost::shared_array<char> what_buffer;

      public:

    stop_info(
      const char* source_file_,
      size_t source_line_,
      ssize_t digits_) throw()
    :
      source_file(source_file_),
      source_line(source_line_),
      digits(digits_),
      what_buffer(0)
    {}

    stop_info(
      const char* source_file_,
      size_t source_line_,
      std::string const& message_) throw()
    :
      source_file(source_file_),
      source_line(source_line_),
      digits(-1),
      message(message_),
      what_buffer(0)
    {}

    ~stop_info() throw() {}

    const char*
    what() const throw()
    {
      if (what_buffer.get() == 0) {
        std::ostringstream o;
        o << "STOP at " << source_file << "(" << source_line << ")";
        if (message.size() != 0) {
          o << ": " << message;
        }
        else if (digits > 0) {
          o << ": code=" << digits;
        }
        std::string result = o.str();
        size_t n = result.size();
        what_buffer.reset(new char[n + 1]);
        std::memcpy(what_buffer.get(), result.data(), n);
        what_buffer[n] = '\0';
      }
      return what_buffer.get();
    }
  };

} // namespace fem

#define FEM_STOP(arg) \
  throw fem::stop_info(__FILE__, __LINE__, arg)

#endif // GUARD
