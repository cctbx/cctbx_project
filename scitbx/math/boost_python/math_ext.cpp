#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/erf.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

namespace scitbx { namespace math {
namespace boost_python { namespace {

  void init_module()
  {
    using namespace boost::python;
    def("erf", (double(*)(double const&)) erf);
    def("erfc", (double(*)(double const&)) erfc);
    def("erfcx", (double(*)(double const&)) erfcx);
  }

}}}} // namespace scitbx::math::boost_python::<anonymous>

BOOST_PYTHON_MODULE(math_ext)
{
  scitbx::math::boost_python::init_module();
}
