/* The main purpose of this module is to provide access to the
   Boost.Python metaclass via meta_ext.empty.__class__.
   See also:
     boost/libs/python/doc/tutorial/doc/quickstart.txt, keyword injector
     boost.python.injector (boost/python.py)
 */

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>

namespace {

  inline
  std::string
  to_str(long value)
  {
    char buf[256];
    sprintf(buf, "%ld", value);
    return std::string(buf);
  }

  std::string
  platform_info()
  {
    std::string result;
#if defined(__alpha__)
    result += " __alpha__";
#endif
#if defined(__host_mips)
    result += " __host_mips";
#endif
#if defined(__i386__)
    result += " __i386__";
#endif
#if defined(__x86_64__)
    result += " __x86_64__";
#endif
#if defined(__linux)
    result += " __linux";
#endif
#if defined(__osf__)
    result += " __osf__";
#endif
#if defined(__sgi)
    result += " __sgi";
#endif
#if defined(_WIN32)
    result += " _WIN32";
#endif
#if defined(_WIN64)
    result += " _WIN64";
#endif
#if defined(__APPLE_CC__)
    result += " __APPLE_CC__=" + to_str(__APPLE_CC__);
#endif
#if defined(_COMPILER_VERSION)
    result += " _COMPILER_VERSION=" + to_str(_COMPILER_VERSION);
#endif
#if defined(__DECCXX_VER)
    result += " __DECCXX_VER=" + to_str(__DECCXX_VER);
#endif
#if defined(__EDG_VERSION__)
    result += " __EDG_VERSION__=" + to_str(__EDG_VERSION__);
#endif
#if defined(__GNUC__)
    result += " __GNUC__=" + to_str(__GNUC__);
#endif
#if defined(__GNUC_MINOR__)
    result += " __GNUC_MINOR__=" + to_str(__GNUC_MINOR__);
#endif
#if defined(__GNUC_PATCHLEVEL__)
    result += " __GNUC_PATCHLEVEL__=" + to_str(__GNUC_PATCHLEVEL__);
#endif
#if defined(__INTEL_COMPILER)
    result += " __INTEL_COMPILER=" + to_str(__INTEL_COMPILER);
#endif
#if defined(_MSC_VER)
    result += " _MSC_VER=" + to_str(_MSC_VER);
#endif
    if (result.size() > 0) result = result.substr(1);
    return result;
  }

} // namespace anonymous

namespace { struct empty {}; }

BOOST_PYTHON_MODULE(boost_python_meta_ext)
{
  using namespace boost::python;
  def("platform_info", platform_info);
  class_<empty>("empty", no_init);
}
