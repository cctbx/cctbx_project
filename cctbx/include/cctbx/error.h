// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
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
#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush
#if 0
using std::cout;
using std::endl;
#endif
#endif

//! Common cctbx namespace.
namespace cctbx {

  //! All cctbx exceptions are derived from this class.
  class error : public std::exception {
    protected:
      std::string message;
    public:
      // NO DEFAULT CONSTRUCTOR
      //! General cctbx error message.
      error(const std::string& msg) throw();
      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, const std::string& msg = "",
            bool Internal = true) throw();
      //! Virtual destructor.
      virtual ~error() throw() { }
      //! Access to the error messages.
      virtual const char* what() const throw();
  };

  //! Special class for "Index out of range." exceptions.
  /*! These exceptions are propagated to Python as IndexError.
   */
  class error_index : public error {
    public:
      //! Default constructor. The message may be customized.
      explicit
      error_index(const std::string& msg = "Index out of range.") throw()
        : error(msg) {}
      //! Virtual destructor.
      virtual ~error_index() throw() {}
  };

} // namespace cctbx

//! For throwing an "Internal Error" exception.
#define cctbx_internal_error() cctbx::error(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define cctbx_not_implemented() cctbx::error(__FILE__, __LINE__, \
             "Not implemented.")
//! Custom cctbx assertion.
#define cctbx_assert(bool) \
  if (!(bool)) throw cctbx::error(__FILE__, __LINE__,\
    "cctbx_assert(" # bool ") failure.")

#endif // CCTBX_ERROR_H
