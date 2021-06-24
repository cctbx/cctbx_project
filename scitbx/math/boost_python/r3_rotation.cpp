#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/math/r3_rotation.h>

namespace scitbx { namespace math { namespace r3_rotation {
namespace {

  struct axis_and_angle_from_matrix_wrappers
  {
    typedef axis_and_angle_from_matrix<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("r3_rotation_axis_and_angle_from_matrix", no_init)
        .def(init<mat3<double> const&>((arg("r"))))
        .add_property("axis", make_getter(&w_t::axis, rbv()))
        .def("angle", &w_t::angle, (arg("deg")=false))
        .def("as_matrix", &w_t::as_matrix)
        .def("as_unit_quaternion", &w_t::as_unit_quaternion)
      ;
    }
  };

}} // namespace <anonymous>::r3_rotation

namespace boost_python {

  void wrap_r3_rotation()
  {
    using namespace boost::python;
    def("r3_rotation_axis_and_angle_as_matrix",
      (mat3<double>(*)(
        vec3<double> const&, double, bool deg, double const&))
      r3_rotation::axis_and_angle_as_matrix, (
        arg("axis"),
        arg("angle"),
        arg("deg")=false,
        arg("min_axis_length")=1e-15));
    def("r3_rotation_axis_and_angle_as_matrix",
      (af::shared<mat3<double> >(*)(
        af::shared<vec3<double> > const&,
        af::shared<double> const&,
        bool,
        double const&)
      ) r3_rotation::axis_and_angle_as_matrix, (
        arg("axes"),
        arg("angle"),
        arg("deg")=false,
        arg("min_axis_length")=1e-15));

    def("r3_rotation_axis_and_angle_as_unit_quaternion",
      (af::tiny<double, 4>(*)(
        vec3<double> const&, double, bool deg, double const&))
      r3_rotation::axis_and_angle_as_unit_quaternion, (
        arg("axis"),
        arg("angle"),
        arg("deg")=false,
        arg("min_axis_length")=1e-15));

    r3_rotation::axis_and_angle_from_matrix_wrappers::wrap();

    def("r3_rotation_vector_to_vector",
      (mat3<double>(*)(
        vec3<double> const&,
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_vector, (
        arg("given_unit_vector"),
        arg("target_unit_vector"),
        arg("sin_angle_is_zero_threshold")=1e-10));

    def("r3_rotation_vector_to_001",
      (mat3<double>(*)(
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_001, (
        arg("given_unit_vector"),
        arg("sin_angle_is_zero_threshold")=1e-10));

    def("r3_rotation_vector_to_010",
      (mat3<double>(*)(
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_010, (
        arg("given_unit_vector"),
        arg("sin_angle_is_zero_threshold")=1e-10));

    def("r3_rotation_vector_to_100",
      (mat3<double>(*)(
        vec3<double> const&,
        double const&))
      r3_rotation::vector_to_100, (
        arg("given_unit_vector"),
        arg("sin_angle_is_zero_threshold")=1e-10));

    def("r3_rotation_unit_quaternion_as_matrix",
      (mat3<double>(*)(af::tiny<double, 4> const&))
        r3_rotation::unit_quaternion_as_matrix, (arg("q")));
    def("r3_rotation_matrix_as_unit_quaternion",
      (af::tiny<double, 4>(*)(mat3<double> const&))
        r3_rotation::matrix_as_unit_quaternion, (arg("r")));
  }

}}} // namespace scitbx::math::boost_python
