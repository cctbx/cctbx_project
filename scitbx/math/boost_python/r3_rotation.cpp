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
        .def(init<mat3<double> const&>((arg_("r"))))
        .add_property("axis", make_getter(&w_t::axis, rbv()))
        .def("angle", &w_t::angle, angle_overloads((arg_("deg")=false)))
        .def("as_matrix", &w_t::as_matrix)
        .def("as_unit_quaternion", &w_t::as_unit_quaternion)
      ;
    }
  };

}} // namespace <anonymous>::r3_rotation

namespace boost_python {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    r3_rotation_axis_and_angle_as_matrix_overloads,
    r3_rotation::axis_and_angle_as_matrix, 2, 4)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    r3_rotation_vector_to_vector_overloads,
    r3_rotation::vector_to_vector, 2, 3)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    r3_rotation_vector_to_001_overloads,
    r3_rotation::vector_to_001, 1, 2)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    r3_rotation_vector_to_010_overloads,
    r3_rotation::vector_to_010, 1, 2)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    r3_rotation_vector_to_100_overloads,
    r3_rotation::vector_to_100, 1, 2)

  void wrap_r3_rotation()
  {
    using namespace boost::python;
    def("r3_rotation_axis_and_angle_as_matrix",
      (mat3<double>(*)(
        vec3<double> const&, double, bool deg, double const&))
      r3_rotation::axis_and_angle_as_matrix,
        r3_rotation_axis_and_angle_as_matrix_overloads((
          arg_("axis"),
          arg_("angle"),
          arg_("deg")=false,
          arg_("min_axis_length")=1.e-15)));

    r3_rotation::axis_and_angle_from_matrix_wrappers::wrap();

    def("r3_rotation_vector_to_vector",
      (mat3<double>(*)(
        vec3<double> const&,
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_vector,
        r3_rotation_vector_to_vector_overloads((
          arg_("given_unit_vector"),
          arg_("target_unit_vector"),
          arg_("sin_angle_is_zero_threshold")=1.e-10)));

    def("r3_rotation_vector_to_001",
      (mat3<double>(*)(
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_001,
        r3_rotation_vector_to_001_overloads((
          arg_("given_unit_vector"),
          arg_("sin_angle_is_zero_threshold")=1.e-10)));

    def("r3_rotation_vector_to_010",
      (mat3<double>(*)(
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_010,
        r3_rotation_vector_to_010_overloads((
          arg_("given_unit_vector"),
          arg_("sin_angle_is_zero_threshold")=1.e-10)));

    def("r3_rotation_vector_to_100",
      (mat3<double>(*)(
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_100,
        r3_rotation_vector_to_100_overloads((
          arg_("given_unit_vector"),
          arg_("sin_angle_is_zero_threshold")=1.e-10)));

    def("r3_rotation_unit_quaternion_as_matrix",
      (mat3<double>(*)(
        double const&, double const&, double const&, double const&))
      r3_rotation::unit_quaternion_as_matrix, ((
        arg_("q0"), arg_("q1"), arg_("q2"), arg_("q3"))));
  }

}}} // namespace scitbx::math::boost_python
