#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>

namespace cctbx { namespace mintbx { namespace boost_python {

  void wrap_k_b_scaling();

namespace {

  void init_module()
  {
    wrap_k_b_scaling();
  }

} // namespace <anonymous>
}}} // namespace cctbx::mintbx::boost_python

BOOST_PYTHON_MODULE(cctbx_mintbx_ext)
{
  cctbx::mintbx::boost_python::init_module();
}
