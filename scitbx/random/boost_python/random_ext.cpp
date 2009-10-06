#include <boost/python/module.hpp>

namespace scitbx { namespace random { namespace boost_python {

  void wrap_random();

  namespace {
    void init_module() {
      wrap_random();
    }
  }
}}}

BOOST_PYTHON_MODULE(scitbx_random_ext)
{
  scitbx::random::boost_python::init_module();
}
