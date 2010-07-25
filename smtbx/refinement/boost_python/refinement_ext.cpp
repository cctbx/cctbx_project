#include <boost/python/module.hpp>

namespace smtbx { namespace refinement {

namespace boost_python {
  void wrap_minimization();

  namespace {
    void init_module() {
      wrap_minimization();
    }
  }

}}} // end namespace smtbx::refinement::boost_python

BOOST_PYTHON_MODULE(smtbx_refinement_ext)
{
  smtbx::refinement::boost_python::init_module();
}
