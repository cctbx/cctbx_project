#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <boost/python/module.hpp>
#include <vector>
#include <set>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  void init_module()
  {
    shared_wrapper<std::vector<std::size_t> >::wrap(
      "std_vector_size_t");
    scitbx::boost_python::container_conversions
    ::tuple_mapping_variable_capacity<std::vector<std::size_t> >();

    shared_wrapper<std::set<std::size_t> >::wrap(
      "std_set_size_t");
    scitbx::boost_python::container_conversions
    ::tuple_mapping_set<std::set<std::size_t> >();
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_array_family_shared_ext)
{
  scitbx::af::boost_python::init_module();
}
