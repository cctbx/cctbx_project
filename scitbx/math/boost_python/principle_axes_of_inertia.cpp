#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <scitbx/math/principle_axes_of_inertia.h>

namespace scitbx { namespace math {
namespace {

  struct principle_axes_of_inertia_wrappers
  {
    typedef principle_axes_of_inertia<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("principle_axes_of_inertia", no_init)
        .def(init<af::const_ref<vec3<double> > const&>((arg_("points"))))
        .def(init<af::const_ref<vec3<double> > const&,
                  af::const_ref<double> const&>(
          (arg_("points"), arg_("weights"))))
        .def("center_of_mass", &w_t::center_of_mass, ccr())
        .def("inertia_tensor", &w_t::inertia_tensor, ccr())
        .def("eigensystem", &w_t::eigensystem, rir())
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_principle_axes_of_inertia()
  {
    principle_axes_of_inertia_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
