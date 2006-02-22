#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/math/r3_rotation.h>

namespace scitbx { namespace math { namespace r3_rotation {
namespace {

  struct axis_and_angle_from_matrix_wrappers
  {
    typedef axis_and_angle_from_matrix<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      angle_overloads, angle, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("r3_rotation_axis_and_angle_from_matrix", no_init)
        .def(init<scitbx::mat3<double> const&>((arg_("r"))))
        .add_property("axis", make_getter(&w_t::axis, rbv()))
        .def("angle", &w_t::angle, angle_overloads((arg_("deg")=false)))
        .def("as_matrix", &w_t::as_matrix)
      ;
    }
  };

}} // namespace <anonymous>::r3_rotation

namespace boost_python {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    r3_rotation_axis_and_angle_as_matrix_overloads,
    r3_rotation::axis_and_angle_as_matrix, 2, 4)

  void wrap_r3_rotation()
  {
    using namespace boost::python;
    def("r3_rotation_axis_and_angle_as_matrix",
      (scitbx::mat3<double>(*)(
        scitbx::vec3<double> const&, double, bool deg, double const&))
      r3_rotation::axis_and_angle_as_matrix,
        r3_rotation_axis_and_angle_as_matrix_overloads((
          arg_("axis"),
          arg_("angle"),
          arg_("deg")=false,
          arg_("min_axis_length")=1.e-15)));

    r3_rotation::axis_and_angle_from_matrix_wrappers::wrap();
  }

}}} // namespace scitbx::math::boost_python
