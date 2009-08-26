#include <boost/python/module.hpp>
#include <omptbx/omp_or_stubs.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_icering();
  void wrap_flex_w_spot();
  void wrap_flex_point();
  void wrap_flex_intxy();

namespace {

  void init_module()
  {
    using namespace boost::python;

    wrap_flex_icering();
    wrap_flex_w_spot();
    wrap_flex_point();
    wrap_flex_intxy();
  }

} // namespace <anonymous>
}}} // namespace scitbx::af::boost_python

BOOST_PYTHON_MODULE(spotfinder_array_family_flex_ext)
{
  scitbx::af::boost_python::init_module();
}
