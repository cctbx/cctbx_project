#include <boost/python/module.hpp>
#include <ccp4_errno.h>

namespace iotbx { namespace mtz { namespace boost_python {

  void wrap_object();
  void wrap_crystal();
  void wrap_dataset();
  void wrap_column();

namespace {

  void init_module()
  {
    CCP4::ccp4_liberr_verbosity(0);
    wrap_object();
    wrap_crystal();
    wrap_dataset();
    wrap_column();
  }

}}}} // namespace iotbx::mtz::boost_python::<anonymous>

BOOST_PYTHON_MODULE(iotbx_mtz_wrapper_ext)
{
  iotbx::mtz::boost_python::init_module();
}
