#include <boost/python/module.hpp>

namespace smtbx { namespace refinement { namespace boost_python {

  void wrap_minimization();
  void wrap_u_cart_special_position_constraints();

  namespace {

    void init_module() {
      wrap_minimization();
      wrap_u_cart_special_position_constraints();
    }

  } // namespace anonymous
}}} // namespace smtbx::refinement::boost_python


BOOST_PYTHON_MODULE(smtbx_refinement_ext)
{
  smtbx::refinement::boost_python::init_module();
}
