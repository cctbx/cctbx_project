#include <boost/python/module.hpp>

namespace cctbx { namespace restraints { namespace boost_python {

  void wrap_bond();
  void wrap_repulsion();
  void wrap_angle();
  void wrap_dihedral();
  void wrap_chirality();
  void wrap_planarity();
  void wrap_pair_proxies();
  void wrap_bonded_interactions();

namespace {

  void init_module()
  {
    wrap_bond();
    wrap_repulsion();
    wrap_angle();
    wrap_dihedral();
    wrap_chirality();
    wrap_planarity();
    wrap_pair_proxies();
    wrap_bonded_interactions();
  }

} // namespace <anonymous>
}}} // namespace cctbx::restraints::boost_python

BOOST_PYTHON_MODULE(cctbx_restraints_ext)
{
  cctbx::restraints::boost_python::init_module();
}
