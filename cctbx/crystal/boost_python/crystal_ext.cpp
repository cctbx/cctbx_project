#include <boost/python/module.hpp>

namespace cctbx { namespace crystal { namespace boost_python {

  void wrap_direct_space_asu();

namespace {

  void init_module()
  {
    wrap_direct_space_asu();
  }

} // namespace <anonymous>

}}} // namespace cctbx::sgtbx::boost_python

BOOST_PYTHON_MODULE(cctbx_crystal_ext)
{
  cctbx::crystal::boost_python::init_module();
}
