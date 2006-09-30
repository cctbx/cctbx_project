/*! \file
    Declarations and macros for exception handling.
 */

#ifndef MMTBX_ERROR_H
#define MMTBX_ERROR_H

#include <stdio.h>
#include <exception>
#include <string>

#define MMTBX_CHECK_POINT std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

//! Common mmtbx namespace.
namespace mmtbx {

  namespace detail {

    inline
    std::string
    compose_error_message(
      std::string const& prefix,
      std::string const& msg)
    {
      return prefix + " Error: " + msg;
    }

    inline
    std::string
    compose_error_message(
      std::string const& prefix,
      const char* file,
      long line,
      std::string const& msg,
      bool internal)
    {
      const char *s = "";
      if (internal) s = " Internal";
      char buf[64];
      sprintf(buf, "%ld", line);
      std::string result = prefix + s + " Error: " + file + "(" + buf + ")";
      if (msg.size()) result += std::string(": ") + msg;
      return result;
    }

  }

  //! All mmtbx exceptions are derived from this class.
  class error : public std::exception
  {
    public:
      //! General mmtbx error message.
      explicit
      error(std::string const& msg) throw()
      : msg_(detail::compose_error_message("mmtbx", msg))
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
      : msg_(detail::compose_error_message(
          "mmtbx", file, line, msg, internal))
      {}

      //! Virtual destructor.
      virtual ~error() throw() {}

      //! Access to the error messages.
      virtual const char* what() const throw()
      {
        return msg_.c_str();
      }

    protected:
      struct complete_msg_tag {};

      error(std::string const& msg, complete_msg_tag)
      : msg_(msg)
      {}

      std::string msg_;
  };

  //! Special class for "Index out of range." exceptions.
  /*! These exceptions are propagated to Python as IndexError.
   */
  class error_index : public error
  {
    public:
      //! Default constructor. The message may be customized.
      explicit
      error_index(std::string const& msg = "Index out of range.") throw()
      : error(msg)
      {}

      //! Virtual destructor.
      virtual ~error_index() throw() {}
  };

} // namespace mmtbx

//! For throwing an "Internal Error" exception.
#define MMTBX_INTERNAL_ERROR() ::mmtbx::error(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define MMTBX_NOT_IMPLEMENTED() ::mmtbx::error(__FILE__, __LINE__, \
             "Not implemented.")
//! Custom mmtbx assertion.
#define MMTBX_ASSERT(bool) \
  if (!(bool)) throw ::mmtbx::error(__FILE__, __LINE__,\
    "MMTBX_ASSERT(" # bool ") failure.")

#endif // MMTBX_ERROR_H
