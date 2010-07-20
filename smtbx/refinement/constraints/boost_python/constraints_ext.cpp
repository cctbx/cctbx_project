#include <boost/python/module.hpp>

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

namespace boost_python {
  void wrap_reparametrisation();
  void wrap_geometrical_hydrogens();

  namespace {
    void init_module() {
      wrap_reparametrisation();
      wrap_geometrical_hydrogens();
    }
  }

}}}} // end namespace smtbx::refinement::constraints::boost_python

BOOST_PYTHON_MODULE(smtbx_refinement_constraints_ext)
{
  smtbx::refinement::constraints::boost_python::init_module();
}
