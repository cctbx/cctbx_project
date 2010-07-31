#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <boost/python.hpp>

#include <scitbx/r3_utils.hpp>

namespace scitbx { namespace r3_utils {
namespace boost_python {

  void init_module()
  {
    using namespace boost::python;
    {
      typedef clash_detector_simple wt;
      class_<wt, boost::noncopyable>("clash_detector_simple", no_init)
        .def(init<unsigned, double>((
          arg("n_sites"),
          arg("threshold"))))
        .def_readonly("threshold_sq", &wt::threshold_sq)
        .def("add_exclusion", &wt::add_exclusion, (arg("i"), arg("j")))
        .def("has_clash", &wt::has_clash, (arg("sites_cart")))
      ;
    }
  }

}}} // namespace scitbx::r3_utils::boost_python

BOOST_PYTHON_MODULE(scitbx_r3_utils_ext)
{
  scitbx::r3_utils::boost_python::init_module();
}
