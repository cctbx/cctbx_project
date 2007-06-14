/*! \file
    Declarations and macros for exception handling.
 */

#ifndef CCTBX_ERROR_H
#define CCTBX_ERROR_H

#include <scitbx/smart_error.h>

//! Common cctbx namespace.
namespace cctbx {

  //! All cctbx exceptions are derived from this class.
  class error : public scitbx::smart_error<error>
  {
    public:

      //! General cctbx error message.
      explicit
      error(std::string const& msg) throw()
        : scitbx::smart_error<error>("cctbx", msg)
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
        : scitbx::smart_error<error>("cctbx", file, line, msg, internal)
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

} // namespace cctbx

//! For throwing an error exception with file name, line number, and message.
#define CCTBX_ERROR(msg) REPORT_ERROR(cctbx::error, msg)
//! For throwing an "Internal Error" exception.
#define CCTBX_INTERNAL_ERROR() REPORT_INTERNAL_ERROR(cctbx::error)
//! For throwing a "Not implemented" exception.
#define CCTBX_NOT_IMPLEMENTED() REPORT_NOT_IMPLEMENTED(cctbx::error)

//! Custom cctbx assertion.
#define CCTBX_ASSERT(assertion)\
  SMART_ASSERT(cctbx::error, CCTBX_ASSERT, assertion)

#endif // CCTBX_ERROR_H
