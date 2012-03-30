#ifndef SCITBX_ERROR_UTILS_H
#define SCITBX_ERROR_UTILS_H

#include <sstream>
#include <string>
#include <exception>
#include <stdexcept>

#ifdef _MSC_VER
# pragma warning(disable:4355)
#endif

namespace scitbx {

  //! Exception which can plug into the SCITBX_ERROR_UTILS_ASSERT trickery
  /*!
    This class is aimed at being inherited from:
    class error : public std::exception, protected scitbx::error_base<error>
    {
      // define constructors
    }
  */
  template<class DerivedError>
  class error_base : public std::exception
  {
    public:
      DerivedError &SCITBX_ERROR_UTILS_ASSERT_A, &SCITBX_ERROR_UTILS_ASSERT_B;

      DerivedError& derived_this() throw();

      error_base(std::string const& prefix, std::string const& msg) throw();

      error_base(
        std::string const& prefix,
        const char* file, long line,
        std::string const& msg = "",
        bool internal = true) throw();

      error_base(error_base const& e) throw();

      virtual ~error_base() throw();

      virtual const char*
      what() const throw();

      template<typename T>
      DerivedError&
      with_current_value(const T& value, const char* label);

    protected:
      std::string msg_;
  };

  template<class DerivedError>
  DerivedError& error_base<DerivedError>
  ::derived_this() throw()
  {
    return static_cast<DerivedError&>(*this);
  }

  template<class DerivedError>
  error_base<DerivedError>
  ::error_base(
    std::string const& prefix,
    std::string const& msg) throw()
  :
    SCITBX_ERROR_UTILS_ASSERT_A(derived_this()),
    SCITBX_ERROR_UTILS_ASSERT_B(derived_this())
  {
    std::ostringstream o;
    o << prefix << " Error: " << msg;
    msg_ = o.str();
  }

  template<class DerivedError>
  error_base<DerivedError>
  ::error_base(
    std::string const& prefix,
    const char* file, long line,
    std::string const& msg,
    bool internal) throw()
  :
    SCITBX_ERROR_UTILS_ASSERT_A(derived_this()),
    SCITBX_ERROR_UTILS_ASSERT_B(derived_this())
  {
    std::ostringstream o;
    o << prefix;
    if (internal) o << " Internal";
    o << " Error: " << file << "(" << line << ")";
    if (msg.size()) o << ": " << msg;
    msg_ = o.str();
  }

  template<class DerivedError>
  error_base<DerivedError>
  ::error_base(error_base const& e) throw()
  :
    std::exception(e),
    SCITBX_ERROR_UTILS_ASSERT_A(derived_this()),
    SCITBX_ERROR_UTILS_ASSERT_B(derived_this())
  {
    msg_ += e.msg_;
  }

  template<class DerivedError>
  error_base<DerivedError>
  ::~error_base() throw() {}

  template<class DerivedError>
  const char*
  error_base<DerivedError>
  ::what() const throw() { return msg_.c_str(); }

  template<class DerivedError>
  template<typename T>
  DerivedError& error_base<DerivedError>
  ::with_current_value(
    const T& value,
    const char* label)
  {
    std::ostringstream o;
    o << "\n" << "  " << label << " = " << value;
    msg_ += o.str();
    return derived_this();
  }

} // namespace scitbx

//! For throwing an error exception with file name, line number, and message.
#define SCITBX_ERROR_UTILS_REPORT(error_class, msg) \
  error_class(__FILE__, __LINE__, msg, false)
//! For throwing an "Internal Error" exception.
#define SCITBX_ERROR_UTILS_REPORT_INTERNAL(error_class) \
  error_class(__FILE__, __LINE__)
//! For throwing a "Not implemented" exception.
#define SCITBX_ERROR_UTILS_REPORT_NOT_IMPLEMENTED(error_class) \
  error_class(__FILE__, __LINE__, "Not implemented.")

//! Custom assertion.
/*!
Here is an example of use:
 \code
 SCITBX_ERROR_UTILS_ASSERT(error_class, n > 0 && m < n)(m)(n);
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

#define SCITBX_ERROR_UTILS_ASSERT_A(x) SCITBX_ERROR_UTILS_ASSERT_OP(x, B)
#define SCITBX_ERROR_UTILS_ASSERT_B(x) SCITBX_ERROR_UTILS_ASSERT_OP(x, A)
#define SCITBX_ERROR_UTILS_ASSERT_OP(x, next) \
  SCITBX_ERROR_UTILS_ASSERT_A.with_current_value( \
    (x), #x).SCITBX_ERROR_UTILS_ASSERT_ ## next

#define SCITBX_ERROR_UTILS_ASSERT(error_class, assert_macro, assertion) \
  if (!(assertion)) throw error_class(__FILE__, __LINE__, \
    #assert_macro"(" # assertion ") failure.").SCITBX_ERROR_UTILS_ASSERT_A

#endif // SCITBX_ERROR_UTILS_H
