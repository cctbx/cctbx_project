#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <scitbx/math/principal_axes_of_inertia.h>

namespace scitbx { namespace math {
namespace {

  struct principal_axes_of_inertia_wrappers
  {
    typedef principal_axes_of_inertia<> w_t;
    typedef principal_axes_of_inertia_2d<> w_t_2d;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("principal_axes_of_inertia", no_init)
        .def(init<af::const_ref<vec3<double> > const&>((arg_("points"))))
        .def(init<af::const_ref<vec3<double> > const&,
                  af::const_ref<double> const&>(
          (arg_("points"), arg_("weights"))))
        .def("center_of_mass", &w_t::center_of_mass, ccr())
        .def("inertia_tensor", &w_t::inertia_tensor, ccr())
        .def("eigensystem", &w_t::eigensystem, rir())
        .def("distance_to_inertia_ellipsoid_surface",
          &w_t::distance_to_inertia_ellipsoid_surface, (
            arg_("unit_direction")))
      ;
      class_<w_t_2d>("principal_axes_of_inertia_2d", no_init)
        .def(init<af::const_ref<vec2<double> > const&>((arg_("points"))))
        .def(init<af::const_ref<vec2<double> > const&,
                  af::const_ref<double> const&>(
          (arg_("points"), arg_("weights"))))
        .def("center_of_mass", &w_t_2d::center_of_mass, ccr())
        .def("inertia_tensor", &w_t_2d::inertia_tensor, ccr())
        .def("eigensystem", &w_t_2d::eigensystem, rir())
        .def("distance_to_inertia_ellipsoid_surface",
          &w_t_2d::distance_to_inertia_ellipsoid_surface, (
            arg_("unit_direction")))
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_principal_axes_of_inertia()
  {
    principal_axes_of_inertia_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
