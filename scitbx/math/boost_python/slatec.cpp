#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/slatec/lib_cpp.h>

#if defined(__linux__) && defined(__GNUC__)
#define SCITBX_MATH_BOOST_PYTHON_CMATH_LGAMMA
#endif

namespace scitbx { namespace math {

namespace {

#if defined(SCITBX_MATH_BOOST_PYTHON_CMATH_LGAMMA)
  inline double
  cmath_lgamma(double const& x) { return lgamma(x); }
#endif

} // namespace <anonymous>

namespace boost_python {

  void wrap_slatec()
  {
    using namespace boost::python;
    slatec_clear_error();
    def("slatec_dgamma", slatec::dgamma, (arg_("x")));
    def("slatec_dlngam", slatec::dlngam, (arg_("x")));
#if defined(SCITBX_MATH_BOOST_PYTHON_CMATH_LGAMMA)
    def("cmath_lgamma", cmath_lgamma, (arg_("x")));
#endif
  }

}}} // namespace scitbx::math::boost_python
