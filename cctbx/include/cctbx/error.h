/*! \file
    Declarations and macros for exception handling.
 */

#ifndef CCTBX_ERROR_H
#define CCTBX_ERROR_H

#include <exception>
#include <string>

//! Common cctbx namespace.
namespace cctbx {

  //! All cctbx exceptions are derived from this class.
  class error : public std::exception
  {
    public:
      //! General cctbx error message.
      explicit
      error(std::string const& msg) throw();

      //! Error message with file name and line number.
      /*! Used by the macros below.
       */
      error(const char* file, long line, std::string const& msg = "",
            bool internal = true) throw();

      //! Virtual destructor.
      virtual ~error() throw();

      //! Access to the error messages.
      virtual const char* what() const throw();

    protected:
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
      error_index(std::string const& msg = "Index out of range.") throw();

      //! Virtual destructor.
      virtual ~error_index() throw();
  };

} // namespace cctbx

//! For throwing an "Internal Error" exception.
#define CCTBX_INTERNAL_ERROR() ::cctbx::error(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define CCTBX_NOT_IMPLEMENTED() ::cctbx::error(__FILE__, __LINE__, \
             "Not implemented.")
//! Custom cctbx assertion.
#define CCTBX_ASSERT(bool) \
  if (!(bool)) throw ::cctbx::error(__FILE__, __LINE__,\
    "CCTBX_ASSERT(" # bool ") failure.")

#endif // CCTBX_ERROR_H
