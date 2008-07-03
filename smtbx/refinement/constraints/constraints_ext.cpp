#include <boost/python/module.hpp>

namespace smtbx { namespace refinement { namespace constraints {

namespace boost_python {
  void wrap_special_positions();
  void wrap_geometric_hydrogen();

  namespace {
    void init_module() {
      wrap_special_positions();
      wrap_geometric_hydrogen();
    }
  }

}}}} // end namespace smtbx::refinement::constraints::boost_python

BOOST_PYTHON_MODULE(smtbx_refinement_constraints_ext)
{
  smtbx::refinement::constraints::boost_python::init_module();
}
