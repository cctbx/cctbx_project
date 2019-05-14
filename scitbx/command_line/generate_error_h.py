from __future__ import absolute_import, division, print_function
from libtbx.utils import write_this_is_auto_generated
import sys, os

def as_string(module_name):
  mn = module_name
  mn_upper = mn.upper()
  return r"""/*! \file
    Declarations and macros for exception handling.
 */

#ifndef %(mn_upper)s_ERROR_H
#define %(mn_upper)s_ERROR_H

#include <scitbx/error_utils.h>

#define %(mn_upper)s_CHECK_POINT \
  std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush
#define %(mn_upper)s_CHECK_POINT_MSG(msg) \
  std::cout << msg << " @ " __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush
#define %(mn_upper)s_EXAMINE(A) \
  std::cout << "variable " << #A << ": " << A << std::endl << std::flush

//! Common %(mn)s namespace.
namespace %(mn)s {

  //! All %(mn)s exceptions are derived from this class.
  class error : public ::scitbx::error_base<error>
  {
    public:

      //! General %(mn)s error message.
      explicit
      error(std::string const& msg) throw()
        : ::scitbx::error_base<error>("%(mn)s", msg)
      {}

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw()
        : ::scitbx::error_base<error>("%(mn)s", file, line, msg, internal)
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

} // namespace %(mn)s

//! For throwing an error exception with file name, line number, and message.
#define %(mn_upper)s_ERROR(msg) \
  SCITBX_ERROR_UTILS_REPORT(%(mn)s::error, msg)
//! For throwing an "Internal Error" exception.
#define %(mn_upper)s_INTERNAL_ERROR() \
  SCITBX_ERROR_UTILS_REPORT_INTERNAL(%(mn)s::error)
//! For throwing a "Not implemented" exception.
#define %(mn_upper)s_NOT_IMPLEMENTED() \
  SCITBX_ERROR_UTILS_REPORT_NOT_IMPLEMENTED(%(mn)s::error)

//! Custom %(mn)s assertion.
#define %(mn_upper)s_ASSERT(assertion) \
  SCITBX_ERROR_UTILS_ASSERT(%(mn)s::error, %(mn_upper)s_ASSERT, assertion)

#endif // %(mn_upper)s_ERROR_H
""" % vars()

def run(args):
  this_command = os.environ["LIBTBX_DISPATCHER_NAME"]
  if (len(args) != 1):
    print("usage: %s module_name > error.h" % this_command)
    return
  write_this_is_auto_generated(f=sys.stdout, file_name_generator=this_command)
  sys.stdout.write(as_string(module_name=args[0]))

if (__name__ == "__main__"):
  run(sys.argv[1:])
