/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy of cctbx/error.h (R.W. Grosse-Kunstleve)
 */

/*! \file
    Declarations and macros for exception handling.
 */

#ifndef SCITBX_ERROR_H
#define SCITBX_ERROR_H

#include <stdio.h>
#include <exception>
#include <string>

#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

//! Common scitbx namespace.
namespace scitbx {

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

  //! All scitbx exceptions are derived from this class.
  class error : public std::exception
  {
    public:
      //! General scitbx error message.
      explicit
      error(std::string const& msg) throw()
      : msg_(detail::compose_error_message("scitbx", msg))
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
      : msg_(detail::compose_error_message(
          "scitbx", file, line, msg, internal))
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

} // namespace scitbx

//! For throwing an "Internal Error" exception.
#define SCITBX_INTERNAL_ERROR() scitbx::error(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define SCITBX_NOT_IMPLEMENTED() scitbx::error(__FILE__, __LINE__, \
             "Not implemented.")
//! Custom scitbx assertion.
#define SCITBX_ASSERT(bool) \
  if (!(bool)) throw scitbx::error(__FILE__, __LINE__,\
    "scitbx_assert(" # bool ") failure.")

#endif // SCITBX_ERROR_H
