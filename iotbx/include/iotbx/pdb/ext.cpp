#include <boost/python/module.hpp>

namespace iotbx { namespace pdb { namespace boost_python {

  void wrap_hierarchy();
  void wrap_input();

namespace {

  void init_module()
  {
    wrap_hierarchy();
    wrap_input();
  }

}}}} // namespace iotbx::pdb::boost_python::<anonymous>

BOOST_PYTHON_MODULE(iotbx_pdb_ext)
{
  iotbx::pdb::boost_python::init_module();
}
