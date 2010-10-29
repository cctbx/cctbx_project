#include <boost/python/module.hpp>

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

namespace boost_python {
  void wrap_reparametrisation();
  void wrap_geometrical_hydrogens();
  void wrap_special_position();
  void wrap_scatterer_parameters();
  void wrap_symmetry_equivalent_site_parameter();
  void wrap_u_eq_dependent_u_iso();

  namespace {
    void init_module() {
      wrap_reparametrisation();
      wrap_geometrical_hydrogens();
      wrap_special_position();
      wrap_scatterer_parameters();
      wrap_symmetry_equivalent_site_parameter();
      wrap_u_eq_dependent_u_iso();
    }
  }

}}}} // end namespace smtbx::refinement::constraints::boost_python

BOOST_PYTHON_MODULE(smtbx_refinement_constraints_ext)
{
  smtbx::refinement::constraints::boost_python::init_module();
}
