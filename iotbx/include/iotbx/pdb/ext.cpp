#include <boost/python/module.hpp>

namespace iotbx { namespace pdb { namespace boost_python {

  void wrap_common_residue_names();
  void wrap_hierarchy();
  void wrap_input();
  void wrap_xray_structure();

namespace {

  void init_module()
  {
    wrap_common_residue_names();
    wrap_hierarchy();
    wrap_input();
    wrap_xray_structure();
  }

}}}} // namespace iotbx::pdb::boost_python::<anonymous>

BOOST_PYTHON_MODULE(iotbx_pdb_ext)
{
  iotbx::pdb::boost_python::init_module();
}
