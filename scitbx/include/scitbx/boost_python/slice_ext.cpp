#include <boost/python/module.hpp>
#include <scitbx/boost_python/slice.h>

namespace scitbx { namespace boost_python {

namespace {

  void init_module()
  {
    using namespace boost::python;
    scitbx::boost_python::slice_from_python();
  }

}}} // namespace scitbx::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_boost_python_slice_ext)
{
  scitbx::boost_python::init_module();
}
