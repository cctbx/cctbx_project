#include <boost/python/module.hpp>

namespace scitbx { namespace lstbx { namespace normal_equations { namespace boost_python {

  void wrap_normal_equations();

  namespace {
    void init_module() {
      wrap_normal_equations();
    }
  }
}}}}

BOOST_PYTHON_MODULE(scitbx_lstbx_normal_equations_ext)
{
  scitbx::lstbx::normal_equations::boost_python::init_module();
}
