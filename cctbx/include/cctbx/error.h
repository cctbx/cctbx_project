/*! \file
    Declarations and macros for exception handling.
 */

#ifndef CCTBX_ERROR_H
#define CCTBX_ERROR_H

#include <sstream>
#include <exception>
#include <string>

//! Common cctbx namespace.
namespace cctbx {

  //! All cctbx exceptions are derived from this class.
  class error : public std::exception
  {
    public:
      error &CCTBX_ASSERT_HACK_A, &CCTBX_ASSERT_HACK_B;

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

      template<typename T>
      error& with_current_value(const T& value, const char* label) {
        std::ostringstream o;
        o << "\t" << label << " = " << value << "\n";
        values_ += o.str();
        return *this;
      }

    protected:
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
      error_index(std::string const& msg = "Index out of range.") throw();

      //! Virtual destructor.
      virtual ~error_index() throw();
  };

} // namespace cctbx

//! For throwing an error exception with file name, line number, and message.
#define CCTBX_ERROR(msg) ::cctbx::error(__FILE__, __LINE__, msg, false)
//! For throwing an "Internal Error" exception.
#define CCTBX_INTERNAL_ERROR() ::cctbx::error(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define CCTBX_NOT_IMPLEMENTED() ::cctbx::error(__FILE__, __LINE__, \
             "Not implemented.")

//! Custom cctbx assertion.
/*!
Here is an example of use:
 \code
 CCTBX_ASSERT(n > 0 && m < n)(m)(n);
 \endcode
 The first parenthesis contains the expression to assert whereas the subsequent
 parentheses contain the variables whose values will be reported upon failure of
 the assertion. The C++ stream system must know how to handle the type of m and
 n through the operator <<.
 If the condition is violated, an instance of cctbx::error is thrown.

 The implementation uses the hacks described in [1].

 [1] A. Alexandrescu and J. Torjo, "Enhancing Assertions",
 C/C++ Users Journal, August 2003
 */

#define CCTBX_ASSERT_HACK_A(x) CCTBX_ASSERT_HACK_OP(x, B)
#define CCTBX_ASSERT_HACK_B(x) CCTBX_ASSERT_HACK_OP(x, A)
#define CCTBX_ASSERT_HACK_OP(x, next) \
CCTBX_ASSERT_HACK_A.with_current_value((x), #x).CCTBX_ASSERT_HACK_ ## next

#define CCTBX_ASSERT(bool) \
if (!(bool)) throw ::cctbx::error(__FILE__, __LINE__,\
                                   "CCTBX_ASSERT(" # bool ") failure.").CCTBX_ASSERT_HACK_A

#endif // CCTBX_ERROR_H
