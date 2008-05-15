/* The main purpose of this module is to provide access to the
   Boost.Python metaclass.
   See also:
     boost/libs/python/doc/tutorial/doc/quickstart.txt, keyword injector
     boost.python.injector (boost/python.py)
 */

#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/cstdint.hpp>

// for number_of_processors()
#if !defined(_MSC_VER)
#include <unistd.h>
#endif
#if defined(__APPLE_CC__)
#include <sys/sysctl.h>
#endif

#if defined(__linux) \
 || defined(__alpha__) \
 || defined(__host_mips) \
 || defined(__APPLE_CC__) \
 || defined(_MSC_VER)
#include <signal.h>
#define BOOST_ADAPTBX_META_EXT_HAVE_SIGNAL_H
#endif

#if defined(__GNUC__)
#include <fenv.h>
#if defined(__linux) \
 || (defined(__APPLE_CC__) && __APPLE_CC__ >= 5465)
#include <execinfo.h>
#define BOOST_ADAPTBX_META_EXT_HAVE_EXECINFO_H
#if defined(__GNUC__) \
 && ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1))) \
 && !defined(__EDG_VERSION__)
#include <cxxabi.h>
#define BOOST_ADAPTBX_META_EXT_HAVE_CXXABI_H
#endif
#endif
#endif

#if defined(__linux)
#include <gnu/libc-version.h>
#define BOOST_ADAPTBX_META_EXT_HAVE_GNU_LIBC-VERSION_H
#endif

namespace {

  long
  number_of_processors()
  {
#if defined(CTL_HW) && defined(HW_NCPU)
    int mib[2];
    mib[0] = CTL_HW;
    mib[1] = HW_NCPU;
    int ncpu;
    size_t len = sizeof(ncpu);
    sysctl(mib, 2, &ncpu, &len, 0, 0);
    return static_cast<long>(ncpu);
#elif defined(_SC_NPROCESSORS_ONLN)
    return sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    return sysconf(_SC_NPROCESSORS_CONF);
#else
    return 0L;
#endif
  }

  bool
  libtbx_introspection_show_stack()
  {
    using namespace boost::python;
    handle<> hdl(allow_null(PyImport_ImportModule("libtbx.introspection")));
    if (!hdl.get()) {
      PyErr_Clear();
      return false;
    }
#if PY_MAJOR_VERSION > 2 || (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION > 4)
    const
#endif
    char* attr_name = "show_stack_true_stderr";
    // test first, just to be maximally fault tolerant
    if (!PyObject_HasAttrString(hdl.get(), attr_name)) {
      return false;
    }
    hdl = handle<>(allow_null(PyObject_GetAttrString(hdl.get(), attr_name)));
    if (!hdl.get()) { // should never be true
      PyErr_Clear();
      return false;
    }
    hdl = handle<>(allow_null(PyObject_CallFunction(hdl.get(), 0)));
    if (!hdl.get()) {
      PyErr_Clear();
      return false;
    }
    return true;
  }

  bool
  boost_adptbx_libc_backtrace(int n_frames_skip=0)
  {
    bool result = false;
#if defined(BOOST_ADAPTBX_META_EXT_HAVE_EXECINFO_H)
    static const int max_frames = 1024;
    void *array[max_frames];
    int size = backtrace(array, max_frames);
    fprintf(stderr, "libc backtrace (%d frames, most recent call last):\n",
      size - n_frames_skip);
    fflush(stderr);
    char **strings = backtrace_symbols(array, size);
    for(int i=size-1;i>=n_frames_skip;i--) {
      char* s = strings[i];
#if defined(BOOST_ADAPTBX_META_EXT_HAVE_CXXABI_H)
      const char* m_bgn = 0;
#if defined(__APPLE_CC__)
      if (strlen(s) >= 52 && strncmp(s+40, "0x", 2) == 0) {
        m_bgn = strchr(s+40, ' ');
      }
#else // __linux
      m_bgn = strchr(s, '(');
#endif
      if (m_bgn != 0) {
        m_bgn++;
        const char* m_end = strchr(m_bgn,
#if defined(__APPLE_CC__)
        ' '
#else // __linux
        '+'
#endif
        );
        long n = m_end - m_bgn;
        if (n > 0) {
          char* mangled = static_cast<char*>(malloc(n+1));
          if (mangled != 0) {
            strncpy(mangled, m_bgn, n);
            mangled[n] = '\0';
            char* demangled = abi::__cxa_demangle(mangled, 0, 0, 0);
            free(mangled);
            if (demangled != 0) {
              long n1 = m_bgn - s;
              long n2 = strlen(demangled);
              long n3 = strlen(m_end);
              char* b = static_cast<char*>(
                malloc(static_cast<size_t>(n1+n2+n3+1)));
              if (b != 0) {
                strncpy(b, s, n1);
                strncpy(b+n1, demangled, n2);
                strncpy(b+n1+n2, m_end, n3);
                b[n1+n2+n3] = '\0';
                s = b;
              }
              free(demangled);
            }
          }
        }
      }
#endif // BOOST_ADAPTBX_META_EXT_HAVE_CXXABI_H
      fprintf(stderr, "  %s\n", s);
      fflush(stderr);
      if (s != strings[i]) free(s);
      result = true;
    }
    free(strings);
#endif // defined(BOOST_ADAPTBX_META_EXT_HAVE_EXECINFO_H)
    return result;
  }

  void
  show_call_stacks_and_exit(const char* what)
  {
    bool have_py_trace = libtbx_introspection_show_stack();
    bool have_libc_trace = boost_adptbx_libc_backtrace(2);
    const char* hint = "sorry, call stacks not available";
    if (have_py_trace && have_libc_trace) {
      hint = "Python and libc call stacks above";
    }
    else if (have_py_trace) {
      hint = "Python call stack above";
    }
    else if (have_libc_trace) {
      hint = "libc call stack above";
    }
    fprintf(stderr, "%s (%s)\n", what, hint);
    fflush(stderr);
    exit(1);
  }

} // namespace anonymous

extern "C" {

  void
  boost_adaptbx_segmentation_fault_backtrace(int)
  {
    show_call_stacks_and_exit("Segmentation fault");
  }

  void
  boost_adaptbx_bus_error_backtrace(int)
  {
    show_call_stacks_and_exit("Bus error");
  }

  void
  boost_adaptbx_floating_point_error_backtrace(int)
  {
    show_call_stacks_and_exit("Floating-point error");
  }

} // extern "C"

namespace {

  inline
  std::string
  to_str(bool value)
  {
    return std::string(value ? "true" : "false");
  }

  inline
  std::string
  to_str(int value)
  {
    char buf[256];
    sprintf(buf, "%d", value);
    return std::string(buf);
  }

  inline
  std::string
  to_str(unsigned int value)
  {
    char buf[256];
    sprintf(buf, "%u", value);
    return std::string(buf);
  }

  inline
  std::string
  to_str(long value)
  {
    char buf[256];
    sprintf(buf, "%ld", value);
    return std::string(buf);
  }

  inline
  std::string
  to_str(unsigned long value)
  {
    char buf[256];
    sprintf(buf, "%lu", value);
    return std::string(buf);
  }

  std::string
  platform_info()
  {
    std::string result;
    std::string nl = "\n";
    result += "__FILE__ = " __FILE__ "\n";
#if defined(__DATE__)
    result += "__DATE__ = " __DATE__ "\n";
#endif
#if defined(__TIME__)
    result += "__TIME__ = " __TIME__ "\n";
#endif
#if defined(__TIMESTAMP__)
    result += "__TIMESTAMP__ = " __TIMESTAMP__ "\n";
#endif
#if defined(__alpha__)
    result += "__alpha__\n";
#endif
#if defined(__host_mips)
    result += "__host_mips\n";
#endif
#if defined(__i386__)
    result += "__i386__\n";
#endif
#if defined(__x86_64__)
    result += "__x86_64__\n";
#endif
#if defined(__linux)
    result += "__linux\n";
#endif
#if defined(__osf__)
    result += "__osf__\n";
#endif
#if defined(__hpux)
    result += "__hpux\n";
#endif
#if defined(__sgi)
    result += "__sgi\n";
#endif
#if defined(_WIN32)
    result += "_WIN32\n";
#endif
#if defined(_WIN64)
    result += "_WIN64\n";
#endif
#if defined(__APPLE_CC__)
    result += "__APPLE_CC__ = " + to_str(__APPLE_CC__) + nl;
#endif
#if defined(_COMPILER_VERSION)
    result += "_COMPILER_VERSION = " + to_str(_COMPILER_VERSION) + nl;
#endif
#if defined(__DECCXX_VER)
    result += "__DECCXX_VER = " + to_str(__DECCXX_VER) + nl;
#endif
#if defined(__HP_aCC)
    result += "__HP_aCC = " + to_str(__HP_aCC) + nl;
#endif
#if defined(__EDG__)
    result += "__EDG__\n";
#endif
#if defined(__EDG_VERSION__)
    result += "__EDG_VERSION__ = " + to_str(__EDG_VERSION__) + nl;
#endif
#if defined(__EDG_ABI_COMPATIBILITY_VERSION)
    result += "__EDG_ABI_COMPATIBILITY_VERSION = "
            + to_str(__EDG_ABI_COMPATIBILITY_VERSION) + nl;
#endif
#if defined(__EDG_IMPLICIT_USING_STD)
    result += "__EDG_IMPLICIT_USING_STD\n";
#endif
#if defined(__EDG_RUNTIME_USES_NAMESPACES)
    result += "__EDG_RUNTIME_USES_NAMESPACES\n";
#endif
#if defined(__GNUC__)
    result += "__GNUC__ = " + to_str(__GNUC__) + nl;
#endif
#if defined(__GNUC_MINOR__)
    result += "__GNUC_MINOR__ = " + to_str(__GNUC_MINOR__) + nl;
#endif
#if defined(__GNUC_PATCHLEVEL__)
    result += "__GNUC_PATCHLEVEL__ = " + to_str(__GNUC_PATCHLEVEL__) + nl;
#endif
#if defined(BOOST_PYTHON_HAVE_CXXABI_CXA_DEMANGLE_IS_BROKEN)
    result += "boost::python::cxxabi_cxa_demangle_is_broken(): "
            + to_str(boost::python::cxxabi_cxa_demangle_is_broken()) + nl;
#endif
#if defined(__GXX_WEAK__)
    result += "__GXX_WEAK__ = " + to_str(__GXX_WEAK__) + nl;
#endif
#if defined(__IEEE_FLOAT)
    result += "__IEEE_FLOAT = " + to_str(__IEEE_FLOAT) + nl;
#endif
#if defined(__INTEL_COMPILER)
    result += "__INTEL_COMPILER = " + to_str(__INTEL_COMPILER) + nl;
#endif
#if defined(__INTEL_COMPILER_BUILD_DATE)
    result += "__INTEL_COMPILER_BUILD_DATE = "
            + to_str(__INTEL_COMPILER_BUILD_DATE) + nl;
#endif
#if defined(__ICC)
    result += "__ICC = " + to_str(__ICC) + nl;
#endif
#if defined(__LP64__)
    result += "__LP64__ = " + to_str(__LP64__) + nl;
#endif
#if defined(_M_IX86)
    result += "_M_IX86 = " + to_str(_M_IX86) + nl;
#endif
#if defined(_MSC_EXTENSIONS)
    result += "_MSC_EXTENSIONS = " + to_str(_MSC_EXTENSIONS) + nl;
#endif
#if defined(_MSC_VER)
    result += "_MSC_VER = " + to_str(_MSC_VER) + nl;
#endif
#if defined(__VERSION__)
    result += "__VERSION__ = " __VERSION__ "\n";
#endif
#if defined(FE_INEXACT)
    result += "FE_INEXACT = " + to_str(FE_INEXACT) + nl;
#endif
#if defined(FE_DIVBYZERO)
    result += "FE_DIVBYZERO = " + to_str(FE_DIVBYZERO) + nl;
#endif
#if defined(FE_UNDERFLOW)
    result += "FE_UNDERFLOW = " + to_str(FE_UNDERFLOW) + nl;
#endif
#if defined(FE_OVERFLOW)
    result += "FE_OVERFLOW = " + to_str(FE_OVERFLOW) + nl;
#endif
#if defined(FE_INVALID)
    result += "FE_INVALID = " + to_str(FE_INVALID) + nl;
#endif
#if defined(FE_ALL_EXCEPT)
    result += "FE_ALL_EXCEPT = " + to_str(FE_ALL_EXCEPT) + nl;
#endif
#if defined(PY_VERSION)
    result += "PY_VERSION = " PY_VERSION "\n";
#endif
#if defined(PYTHON_API_VERSION)
    result += "PYTHON_API_VERSION = " + to_str(PYTHON_API_VERSION) + nl;
#endif
#undef P
#define P(T) result += "sizeof(" #T ") = " + to_str(sizeof(T)) + nl;
    P(bool)
    P(short)
    P(int)
    P(long)
    P(std::size_t)
    P(void*)
#if !defined(_MSC_VER) || _MSC_VER > 1200
    P(long long)
#endif
    P(float)
    P(double)
    P(long double)
    P(boost::int32_t)
    P(boost::uint32_t)
#if defined(HAVE_WCHAR_H)
    P(wchar_t)
#endif
#if defined(PY_UNICODE_TYPE)
    P(PY_UNICODE_TYPE)
#endif
#undef P
#if defined(BOOST_ADAPTBX_META_EXT_HAVE_GNU_LIBC)
    result += "gnu libc version: ";
    result += gnu_get_libc_version() + nl;
#endif
    return result;
  }

  std::size_t
  sizeof_void_ptr() { return sizeof(void*); }

  void
  enable_signals_backtrace_if_possible()
  {
#if defined(BOOST_ADAPTBX_META_EXT_HAVE_SIGNAL_H)
#if defined(SIGSEGV)
    signal(SIGSEGV, boost_adaptbx_segmentation_fault_backtrace);
#endif
#if defined(SIGBUS)
    signal(SIGBUS, boost_adaptbx_bus_error_backtrace);
#endif
#if defined(SIGFPE)
    signal(SIGFPE, boost_adaptbx_floating_point_error_backtrace);
#endif
#endif
  }

  void
  enable_floating_point_exceptions_if_possible(
    bool
#if defined(FE_DIVBYZERO)
    divbyzero
#endif
    ,
    bool
#if defined(FE_INVALID)
    invalid
#endif
    ,
    bool
#if defined(FE_OVERFLOW)
    overflow
#endif
    )
  {
    int flags = 0;
#if defined(FE_DIVBYZERO)
    if (divbyzero) flags |= FE_DIVBYZERO;
#endif
#if defined(FE_INVALID)
    if (invalid) flags |= FE_INVALID;
#endif
#if defined(FE_OVERFLOW)
    if (overflow) flags |= FE_OVERFLOW;
#endif
    if (flags != 0) {
#if defined(__linux)
      feenableexcept(flags);
#endif
    }
  }

  char
  dereference_char_pointer(const char* pointer) { return *pointer; }

  double
  divide_doubles(double const& x, double const& y) { return x / y; }

} // namespace anonymous

namespace boost_python_meta_ext { struct holder {}; }

BOOST_PYTHON_MODULE(boost_python_meta_ext)
{
  using namespace boost::python;
  def("number_of_processors", number_of_processors);
  def("boost_adptbx_libc_backtrace", boost_adptbx_libc_backtrace);
  def("libtbx_introspection_show_stack", libtbx_introspection_show_stack);
  def("platform_info", platform_info);
  def("sizeof_void_ptr", sizeof_void_ptr);
  def("enable_signals_backtrace_if_possible",
       enable_signals_backtrace_if_possible);
  def("enable_floating_point_exceptions_if_possible",
       enable_floating_point_exceptions_if_possible, (
    arg_("divbyzero"),
    arg_("invalid"),
    arg_("overflow")));
  def("dereference_char_pointer", dereference_char_pointer);
  def("divide_doubles", divide_doubles);
  class_<boost_python_meta_ext::holder>("holder").enable_pickling();
  class_<docstring_options, boost::noncopyable>("docstring_options", no_init)
    .def(init<bool, bool>((
      arg_("show_user_defined"),
      arg_("show_signatures"))))
  ;
}
