#include <boost/python/module.hpp>


namespace cctbx { namespace adp_restraints { namespace boost_python {

  void wrap_isotropic_adp();
  void wrap_rigid_bond();
  void wrap_adp_similarity();
  void wrap_aniso_restraints();
  void wrap_fixed_u_eq_adp();
  void wrap_adp_restraint_base();

namespace {

  void init_module() {
    // order matters: the base must be first
    wrap_adp_restraint_base();
    wrap_isotropic_adp();
    wrap_rigid_bond();
    wrap_adp_similarity();
    wrap_aniso_restraints();
    wrap_fixed_u_eq_adp();
  }

} // namespace <anonymous>
}}} // namespace cctbx::adp_restraints::boost_python

BOOST_PYTHON_MODULE(cctbx_adp_restraints_ext) {
  cctbx::adp_restraints::boost_python::init_module();
}
