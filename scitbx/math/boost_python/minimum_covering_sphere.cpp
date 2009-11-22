#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <scitbx/math/minimum_covering_sphere.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(
  scitbx::math::sphere_3d<>)
SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(
  scitbx::math::minimum_covering_sphere_3d<>)

namespace scitbx { namespace math {
namespace {

  struct sphere_3d_wrappers
  {
    typedef sphere_3d<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("sphere_3d", no_init)
        .def(init<vec3<double> const&, double>(
          (arg("center"), arg("radius"))))
        .def("center", &w_t::center)
        .def("radius", &w_t::radius)
        .def("expand", &w_t::expand, (arg("additional_radius")))
        .def("expand_relative", &w_t::expand_relative,
          (arg("additional_relative_radius")))
        .def("is_inside", &w_t::is_inside, (arg("point")))
        .def("box_min", &w_t::box_min)
        .def("box_max", &w_t::box_max)
      ;
    }
  };

  struct minimum_covering_sphere_3d_wrappers
  {
    typedef minimum_covering_sphere_3d<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<sphere_3d<> > >("minimum_covering_sphere_3d", no_init)
        .def(init<
          af::const_ref<vec3<double> > const&,
          optional<
            double const&,
            double const&,
            vec3<double> const&> >((
              arg("points"),
              arg("epsilon")=1.e-6,
              arg("radius_if_one_or_no_points")=1,
              arg("center_if_no_points")=w_t::default_center_if_no_points)))
        .def("n_iterations", &w_t::n_iterations)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_minimum_covering_sphere()
  {
    sphere_3d_wrappers::wrap();
    minimum_covering_sphere_3d_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
