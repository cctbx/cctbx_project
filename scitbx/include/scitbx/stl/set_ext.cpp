#include <scitbx/stl/set_wrapper.h>
#include <boost/python/module.hpp>

namespace scitbx { namespace stl { namespace boost_python {
namespace {

  void init_module()
  {
    set_wrapper<unsigned>::wrap("unsigned");
  }

}}}} // namespace scitbx::stl::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_stl_set_ext)
{
  scitbx::stl::boost_python::init_module();
}
