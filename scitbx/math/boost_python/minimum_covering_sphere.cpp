#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <scitbx/math/minimum_covering_sphere.h>

namespace scitbx { namespace math {
namespace {

  struct sphere_3d_wrappers
  {
    typedef sphere_3d<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("sphere_3d", no_init)
        .def(init<vec3<double> const&, double>(
          (arg_("center"), arg_("radius"))))
        .def("center", &w_t::center)
        .def("radius", &w_t::radius)
        .def("expand", &w_t::expand, (arg_("additional_radius")))
        .def("is_inside", &w_t::is_inside, (arg_("point")))
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
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t, bases<sphere_3d<> > >("minimum_covering_sphere_3d", no_init)
        .def(init<af::const_ref<vec3<double> > const&,
                  optional<double> >((arg_("points"), arg_("epsilon"))))
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
