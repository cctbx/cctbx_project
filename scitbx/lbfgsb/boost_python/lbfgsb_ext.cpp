#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/lbfgsb/raw.h>
#include <boost/python/module.hpp>

namespace scitbx { namespace lbfgsb { namespace {

  void init_module()
  {
    using namespace boost::python;
  }

}}} // namespace scitbx::lbfgs::<anonymous>

BOOST_PYTHON_MODULE(scitbx_lbfgs_ext)
{
  scitbx::lbfgsb::init_module();
}
