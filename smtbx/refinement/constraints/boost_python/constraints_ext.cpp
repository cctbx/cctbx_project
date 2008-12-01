#include <smtbx/refinement/constraints/boost_python/wrappers.h>
#include <boost/python/module.hpp>

#include <cctbx/xray/scatterer.h>
#include <smtbx/refinement/constraints/special_positions.h>

namespace smtbx { namespace refinement { namespace constraints {

namespace boost_python {
  void wrap_geometric_hydrogen();

  namespace {
    void init_module() {
      constrained_scatterers_wrapper<
        special_positions<double, xray::scatterer<>,
                          af::boost_python::flex_1d> >::wrap(
          "special_positions");
      wrap_geometric_hydrogen();
    }
  }

}}}} // end namespace smtbx::refinement::constraints::boost_python

BOOST_PYTHON_MODULE(smtbx_refinement_constraints_ext)
{
  smtbx::refinement::constraints::boost_python::init_module();
}
