#include <boost/python/module.hpp>

namespace cctbx { namespace crystal { namespace boost_python {

  void wrap_close_packing();
  void wrap_direct_space_asu();
  void wrap_neighbors();
  void wrap_pair_tables();
  void wrap_symmetry();

namespace {

  void init_module()
  {
    wrap_close_packing();
    wrap_direct_space_asu();
    wrap_neighbors();
    wrap_pair_tables();
    wrap_symmetry();
  }

} // namespace <anonymous>

}}} // namespace cctbx::sgtbx::boost_python

BOOST_PYTHON_MODULE(cctbx_crystal_ext)
{
  cctbx::crystal::boost_python::init_module();
}
