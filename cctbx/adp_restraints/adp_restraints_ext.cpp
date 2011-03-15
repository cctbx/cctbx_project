#include <boost/python/module.hpp>


namespace cctbx { namespace adp_restraints { namespace boost_python {

  void wrap_isotropic_adp();
  void wrap_rigid_bond();
  void wrap_adp_similarity();
  void wrap_aniso_restraints();

namespace {

  void init_module()
  {
    wrap_isotropic_adp();
    wrap_rigid_bond();
    wrap_adp_similarity();
    wrap_aniso_restraints();
  }

} // namespace <anonymous>
}}} // namespace cctbx::adp_restraints::boost_python

BOOST_PYTHON_MODULE(cctbx_adp_restraints_ext)
{
  cctbx::adp_restraints::boost_python::init_module();
}
