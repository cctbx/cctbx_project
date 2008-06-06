#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <rstbx/dps_core/direction.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_Direction()
  {
    flex_wrapper<labelit::dptbx::Direction>::plain("Direction");
  }
  void wrap_shared_double_array()
  {
    flex_wrapper<scitbx::af::shared<double> >::plain("flex_double");
    flex_wrapper<scitbx::af::shared<scitbx::vec3<double> > >::plain("flex_vec3double");
  }

}}} // namespace scitbx::af::boost_python
