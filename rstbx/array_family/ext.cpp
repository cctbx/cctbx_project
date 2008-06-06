#include <boost/python/module.hpp>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_Direction();
  void wrap_shared_double_array();

namespace {

  void init_module()
  {
    using namespace boost::python;

    wrap_flex_Direction();
    wrap_shared_double_array();
  }

} // namespace <anonymous>
}}} // namespace scitbx::af::boost_python

BOOST_PYTHON_MODULE(rstbx_array_family_flex_ext)
{
  scitbx::af::boost_python::init_module();
}
