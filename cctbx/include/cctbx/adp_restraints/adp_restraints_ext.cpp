#include <boost/python/module.hpp>


namespace cctbx { namespace adp_restraints { namespace boost_python {

  void wrap_rigid_bond();

namespace {

  void init_module()
  {
    wrap_rigid_bond();
  }

} // namespace <anonymous>
}}} // namespace cctbx::adp_restraints::boost_python

BOOST_PYTHON_MODULE(cctbx_adp_restraints_ext)
{
  cctbx::adp_restraints::boost_python::init_module();
}
