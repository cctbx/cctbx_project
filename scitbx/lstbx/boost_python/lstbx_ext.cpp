#include <boost/python/module.hpp>

namespace scitbx { namespace lstbx { namespace boost_python {

  void wrap_normal_equations();

  namespace {
    void init_module() {
      wrap_normal_equations();
    }
  }
}}}

BOOST_PYTHON_MODULE(scitbx_lstbx_ext)
{
  scitbx::lstbx::boost_python::init_module();
}
