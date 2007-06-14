/*! \file
    Declarations and macros for exception handling.
 */

#ifndef MMTBX_ERROR_H
#define MMTBX_ERROR_H

#include <scitbx/smart_error.h>

//! Common mmtbx namespace.
namespace mmtbx {

  //! All mmtbx exceptions are derived from this class.
  class error : public scitbx::smart_error<error>
  {
    public:

      //! General mmtbx error message.
      explicit
      error(std::string const& msg) throw()
        : scitbx::smart_error<error>("mmtbx", msg)
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
        : scitbx::smart_error<error>("mmtbx", file, line, msg, internal)
      {}
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
  };

} // namespace mmtbx

//! For throwing an error exception with file name, line number, and message.
#define MMTBX_ERROR(msg) REPORT_ERROR(mmtbx::error, msg)
//! For throwing an "Internal Error" exception.
#define MMTBX_INTERNAL_ERROR() REPORT_INTERNAL_ERROR(mmtbx::error)
//! For throwing a "Not implemented" exception.
#define MMTBX_NOT_IMPLEMENTED() REPORT_NOT_IMPLEMENTED(mmtbx::error)

//! Custom mmtbx assertion.
#define MMTBX_ASSERT(assertion)\
  SMART_ASSERT(mmtbx::error, MMTBX_ASSERT, assertion)

#endif // MMTBX_ERROR_H
