// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

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
  if (!(bool)) throw cctbx::error(__FILE__, __LINE__,\
    "cctbx_assert(" # bool ") failure.")

#endif // CCTBX_ERROR_H
