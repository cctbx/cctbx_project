#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <boost/python/module.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <vector>
#include <set>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  void init_module()
  {
    typedef boost::python::return_internal_reference<> rir;
    shared_wrapper<std::vector<unsigned>, rir>::wrap("stl_vector_unsigned");
    shared_wrapper<std::set<unsigned>, rir>::wrap("stl_set_unsigned");
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_array_family_shared_ext)
{
  scitbx::af::boost_python::init_module();
}
