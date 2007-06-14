#ifndef SMART_ERROR_FORMAT_H
#define SMART_ERROR_FORMAT_H

#include <sstream>
#include <exception>
#include <string>

#ifdef _MSC_VER
# pragma warning(disable:4355)
#endif

namespace scitbx {

  //! Exception which can plug into the SMART_ASSERT trickery
  /*!
    This class is aimed at being inherited from:
    class error : public std::exception, protected smart_error<error>
    {
      // define constructors
    }
  */
  template<class DerivedError>
  class smart_error : public std::exception
  {
    public:
      DerivedError &SMART_ASSERT_A, &SMART_ASSERT_B;

      DerivedError& derived_this() throw() {
        return static_cast<DerivedError&>(*this);
      }

      smart_error(std::string const& prefix, std::string const& msg) throw()
        : SMART_ASSERT_A(derived_this()), SMART_ASSERT_B(derived_this())
      {
        o << prefix << " Error: " << msg;
      }

      smart_error(std::string const& prefix,
                   const char* file, long line, std::string const& msg = "",
                   bool internal = true) throw()
        : SMART_ASSERT_A(derived_this()), SMART_ASSERT_B(derived_this())
      {
        o << prefix;
        if (internal) o << " Internal";
        o << "Error: " << file << "(" << line << ")";
        if (msg.size()) o << ": " << msg;
      }

      smart_error(smart_error const& e) throw()
        : SMART_ASSERT_A(derived_this()), SMART_ASSERT_B(derived_this())
      {
        o << e.o.str();
      }

      virtual ~smart_error() throw() {}

      virtual const char* what() const throw()
      {
        return o.str().c_str();
      }

      template<typename T>
      DerivedError& with_current_value(const T& value, const char* label)
      {
        o << "\n" << "\t" << label << " = " << value;
        return derived_this();
      }

    protected:
      std::ostringstream o;
  };

}

//! For throwing an error exception with file name, line number, and message.
#define REPORT_ERROR(error_class, msg)\
  error_class(__FILE__, __LINE__, msg, false)
//! For throwing an "Internal Error" exception.
#define REPORT_INTERNAL_ERROR(error_class)\
  error_class(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define REPORT_NOT_IMPLEMENTED(error_class)\
  error_class(__FILE__, __LINE__, "Not implemented.")

//! Custom assertion.
/*!
Here is an example of use:
 \code
 SMART_ASSERT(error_class, n > 0 && m < n)(m)(n);
 \endcode
 The first parenthesis contains the expression to assert whereas the subsequent
 parentheses contain the variables whose values will be reported upon failure of
 the assertion. The C++ stream system must know how to handle the type of m and
 n through the operator <<.
 If the condition is violated, an instance of error_class is thrown. That
 class should inherit from error_format<namespace::error> as explained in the
 documentation of error_format.

 The implementation uses the tricks described in [1].

 [1] A. Alexandrescu and J. Torjo, "Enhancing Assertions",
     C/C++ Users Journal, August 2003
 */

#define SMART_ASSERT_A(x) SMART_ASSERT_OP(x, B)
#define SMART_ASSERT_B(x) SMART_ASSERT_OP(x, A)
#define SMART_ASSERT_OP(x, next) \
SMART_ASSERT_A.with_current_value((x), #x).SMART_ASSERT_ ## next

#define SMART_ASSERT(error_class, assert_macro, assertion) \
if (!(assertion)) throw error_class(__FILE__, __LINE__,\
                      #assert_macro"(" # assertion ") failure.").SMART_ASSERT_A

#endif
