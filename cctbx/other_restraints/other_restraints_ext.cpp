#include <boost/python/module.hpp>


namespace cctbx { namespace other_restraints { namespace boost_python {

  void wrap_sump();

namespace {

  void init_module() {
    wrap_sump();
  }

} // namespace <anonymous>
}}} // namespace cctbx::adp_restraints::boost_python

BOOST_PYTHON_MODULE(cctbx_other_restraints_ext) {
  cctbx::other_restraints::boost_python::init_module();
}
