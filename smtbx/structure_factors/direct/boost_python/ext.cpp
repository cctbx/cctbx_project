#include <boost/python/module.hpp>

namespace smtbx { namespace structure_factors { namespace direct {

  namespace boost_python {

    void wrap_standard_xray();

    namespace {
      void init_module() {
        wrap_standard_xray();
      }
    }
  } // boost_python
}}}

BOOST_PYTHON_MODULE(smtbx_structure_factors_direct_ext)
{
  smtbx::structure_factors::direct::boost_python::init_module();
}
