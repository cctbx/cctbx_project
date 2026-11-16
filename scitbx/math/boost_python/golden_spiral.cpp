#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/vec3.h>
#include <scitbx/math/golden_spiral.h>
#include <boost_adaptbx/boost_range_python.hpp>

namespace scitbx { namespace math {

namespace {

  struct golden_spiral_wrappers
  {
    typedef golden_spiral<scitbx::vec3<double> > w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      boost_adaptbx::python::generic_range_wrapper<w_t::storage_type>::wrap(
        "golden_spiral_points_range");
      class_<w_t>("golden_spiral", no_init)
        .def(init<size_t>((arg("count"))))
        .add_property("unit_area", make_function(&w_t::unit_area,
          return_value_policy<copy_const_reference>()))
        .add_property("count", make_function(&w_t::count,
          return_value_policy<copy_const_reference>()))
        .add_property("points", make_function(&w_t::points,
          return_internal_reference<>()))
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_golden_spiral()
  {
    golden_spiral_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
