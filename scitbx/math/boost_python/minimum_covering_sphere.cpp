#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <scitbx/math/minimum_covering_sphere.h>

namespace scitbx { namespace math {
namespace {

  struct minimum_covering_sphere_wrappers
  {
    typedef minimum_covering_sphere<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("minimum_covering_sphere", no_init)
        .def(init<af::const_ref<vec3<double> > const&,
                  optional<double> >((arg("points"), arg("epsilon"))))
        .def("n_iterations", &w_t::n_iterations)
        .def("center", &w_t::center)
        .def("radius", &w_t::radius)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_minimum_covering_sphere()
  {
    minimum_covering_sphere_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
