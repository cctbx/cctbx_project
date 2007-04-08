/*! \file
    Declarations and macros for exception handling.
 */

#ifndef MMTBX_ERROR_H
#define MMTBX_ERROR_H

#include <scitbx/error.h>

//! Common mmtbx namespace.
namespace mmtbx {

  //! All mmtbx exceptions are derived from this class.
  class error : public std::exception
  {
    public:
      //! General mmtbx error message.
      explicit
      error(std::string const& msg) throw()
      : msg_(scitbx::detail::compose_error_message("mmtbx", msg))
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
      : msg_(scitbx::detail::compose_error_message(
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

} // namespace mmtbx

//! For throwing an error exception with file name, line number, and message.
#define MMTBX_ERROR(msg) ::mmtbx::error(__FILE__, __LINE__, msg, false)
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
