#include <boost/python/module.hpp>

namespace cctbx { namespace dmtbx { namespace boost_python {

  void wrap_triplet_generator();
  void wrap_triplet_phase_relation();

namespace {

  void init_module()
  {
    wrap_triplet_generator();
    wrap_triplet_phase_relation();
  }

} // namespace <anonymous>
}}} // namespace cctbx::dmtbx::boost_python

BOOST_PYTHON_MODULE(cctbx_dmtbx_ext)
{
  cctbx::dmtbx::boost_python::init_module();
}
