/*! \file
    Declarations and macros for exception handling.
 */

#ifndef SCITBX_ERROR_H
#define SCITBX_ERROR_H

#include <stdio.h>
#include <exception>
#include <string>
#include <sstream>

#define SCITBX_CHECK_POINT std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush
#define SCITBX_EXAMINE(A) std::cout << "variable " <<#A<< ": " << A << " " << std::endl << std::flush

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
      std::ostringstream o;
      o << prefix;
      if (internal) o << " Internal";
      o << "Error: " << file << "(" << line << ")";
      if (msg.size()) o << ": " << msg;
      return o.str();
//      const char *s = "";
//      if (internal) s = " Internal";
//      char buf[64];
//      sprintf(buf, "%ld", line);
//      std::string result = prefix + s + " Error: " + file + "(" + buf + ")";
//      if (msg.size()) result += std::string(": ") + msg;
//      return result;
    }

  }

  //! All scitbx exceptions are derived from this class.
  class error : public std::exception
  {
    public:
      error &SCITBX_ASSERT_HACK_A, &SCITBX_ASSERT_HACK_B;

      //! General scitbx error message.
      explicit
      error(std::string const& msg) throw()
      : msg_(detail::compose_error_message("scitbx", msg)),
        SCITBX_ASSERT_HACK_A(*this), SCITBX_ASSERT_HACK_B(*this)
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
      : msg_(detail::compose_error_message(
          "scitbx", file, line, msg, internal)),
        SCITBX_ASSERT_HACK_A(*this), SCITBX_ASSERT_HACK_B(*this)
      {}

      //! Virtual destructor.
      virtual ~error() throw() {}

      //! Access to the error messages.
      virtual const char* what() const throw()
      {
        std::string result = msg_;
        if (values_.length() > 0) result += "\n" + values_;
        return result.c_str();
      }

      template<typename T>
      error& with_current_value(const T& value, const char* label) {
        std::ostringstream o;
        o << "\t" << label << " = " << value << "\n";
        values_ += o.str();
        return *this;
      }

    protected:
      struct complete_msg_tag {};

      error(std::string const& msg, complete_msg_tag)
      : msg_(msg),
        SCITBX_ASSERT_HACK_A(*this), SCITBX_ASSERT_HACK_B(*this)
      {}

      std::string msg_;
      std::string values_;
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

//! For throwing an error exception with file name, line number, and message.
#define SCITBX_ERROR(msg) ::scitbx::error(__FILE__, __LINE__, msg, false)
//! For throwing an "Internal Error" exception.
#define SCITBX_INTERNAL_ERROR() ::scitbx::error(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define SCITBX_NOT_IMPLEMENTED() ::scitbx::error(__FILE__, __LINE__, \
             "Not implemented.")

//! Custom scitbx assertion.
/*!
Here is an example of use:
 \code
 SCITBX_ASSERT(n > 0 && m < n)(m)(n);
 \endcode
 The first parenthesis contains the expression to assert whereas the subsequent
 parentheses contain the variables whose values will be reported upon failure of
 the assertion. The C++ stream system must know how to handle the type of m and
 n through the operator <<.
 If the condition is violated, an instance of scitbx::error is thrown.

 The implementation uses the hacks described in [1].

 [1] A. Alexandrescu and J. Torjo, "Enhancing Assertions",
     C/C++ Users Journal, August 2003
 */

#define SCITBX_ASSERT_HACK_A(x) SCITBX_ASSERT_HACK_OP(x, B)
#define SCITBX_ASSERT_HACK_B(x) SCITBX_ASSERT_HACK_OP(x, A)
#define SCITBX_ASSERT_HACK_OP(x, next) \
SCITBX_ASSERT_HACK_A.with_current_value((x), #x).SCITBX_ASSERT_HACK_ ## next

#define SCITBX_ASSERT(bool) \
  if (!(bool)) throw ::scitbx::error(__FILE__, __LINE__,\
    "SCITBX_ASSERT(" # bool ") failure.").SCITBX_ASSERT_HACK_A

#endif // SCITBX_ERROR_H
