// $Id$

/*! \file
    Declarations and macros for exception handling.
 */

#ifndef CCTBX_ERROR_H
#define CCTBX_ERROR_H

#include <exception>
#include <string>

#if 0
#include <iostream>
using std::cout;
using std::endl;
#define CheckPoint cout << __FILE__ << "(" << __LINE__ << ")" << endl
#endif

//! Common cctbx namespace.
namespace cctbx {

  class error : public std::exception {
    protected:
      std::string message;
    public:
      error(const std::string& msg) throw();
      error(const char* file, long line, const std::string& msg = "") throw();
      virtual ~error() throw() { }
      virtual const char* what() const throw();
  };

} // namespace cctbx

#define cctbx_internal_error() cctbx::error(__FILE__, __LINE__)
#define cctbx_not_implemented() cctbx::error(__FILE__, __LINE__, \
             "Not implemented.")
#define cctbx_assert(bool) \
  if (! (bool)) throw cctbx::error(__FILE__, __LINE__,\
    "cctbx_assert(" # bool ") failure.")

#endif // CCTBX_ERROR_H
