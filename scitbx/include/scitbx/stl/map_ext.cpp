#include <boost/python/module.hpp>
#include <scitbx/stl/map_wrapper.h>
#include <map>

namespace scitbx { namespace stl { namespace boost_python {
namespace {

  void init_module()
  {
    map_wrapper<std::map<std::string, double> >::wrap("stl_string_double");
  }

}}}} // namespace scitbx::stl::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_stl_map_ext)
{
  scitbx::stl::boost_python::init_module();
}
