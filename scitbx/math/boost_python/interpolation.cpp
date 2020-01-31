#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/interpolation.h>
#include <scitbx/math/linear_interpolation.h>
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

  template <typename FloatType>
  void wrap_interp_2d()
  {
    using namespace boost::python;
    def("linear_interpolation_2d",
      (FloatType (*)(
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&,
          FloatType const&)) linear_interpolation_2d, (
            arg("x1"),
            arg("y1"),
            arg("x2"),
            arg("y2"),
            arg("v1"),
            arg("v2"),
            arg("v3"),
            arg("v4"),
            arg("xx"),
            arg("yy")));
  }

} // namespace <anonymous>

namespace boost_python {

  void wrap_interpolation()
  {
    wrap_splines< scitbx::vec2<double> >();
    wrap_splines< scitbx::vec3<double> >();
    wrap_interp_2d <double> ();
  }

}}} // namespace scitbx::math::boost_python
