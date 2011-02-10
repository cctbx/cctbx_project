#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/interpolation.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>

namespace scitbx { namespace math {

namespace {

  template <typename PointType>
  void wrap_splines()
  {
    using namespace boost::python;
    def("interpolate_catmull_rom_spline",
      (af::shared<PointType>(*)(
        PointType const&,
        PointType const&,
        PointType const&,
        PointType const&,
        unsigned)) interpolate_catmull_rom_spline, (
          arg("p0"),
          arg("p1"),
          arg("p2"),
          arg("p3"),
          arg("n_points")));
  }

} // namespace <anonymous>

namespace boost_python {

  void wrap_interpolation()
  {
    wrap_splines< scitbx::vec2<double> >();
    wrap_splines< scitbx::vec3<double> >();
  }

}}} // namespace scitbx::math::boost_python
