#include <boost/python/module.hpp>
#include <scitbx/stl/vector_wrapper.h>
#include <set>

namespace scitbx { namespace stl { namespace boost_python {
namespace {

  void init_module()
  {
    vector_wrapper<unsigned>::wrap("unsigned");

    typedef boost::python::return_internal_reference<> rir;
    vector_wrapper<std::set<unsigned>, rir>::wrap("set_unsigned");
  }

}}}} // namespace scitbx::stl::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_stl_vector_ext)
{
  scitbx::stl::boost_python::init_module();
}
