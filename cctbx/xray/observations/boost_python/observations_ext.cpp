#include <boost/python/module.hpp>

namespace cctbx { namespace xray {
  namespace boost_python {

  void wrap_observations();
  namespace {
    void init_module() {
      wrap_observations();
    }
  }
}}}

BOOST_PYTHON_MODULE(cctbx_xray_observations_ext) {
  cctbx::xray::boost_python::init_module();
}
