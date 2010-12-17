
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/interpolation.h>

namespace scitbx { namespace math {

namespace {

  struct interpolation_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      def("interpolate_catmull_rom_spline", interpolate_catmull_rom_spline, (
        arg("p0"),
        arg("p1"),
        arg("p2"),
        arg("p3")));
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_interpolation()
  {
    interpolation_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
