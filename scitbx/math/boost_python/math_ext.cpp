#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/math/erf.h>
#include <scitbx/math/bessel.h>
#include <scitbx/math/eigensystem.h>
#include <scitbx/math/matrix_inversion.h>
#include <scitbx/math/phase_error.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace math {
namespace boost_python {

  void wrap_gaussian();
  void wrap_golay();
  void wrap_minimum_covering_sphere();
  void wrap_principal_axes_of_inertia();

namespace {

  struct eigensystem_real_symmetric_wrappers
  {
    typedef eigensystem::real_symmetric<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("eigensystem_real_symmetric", no_init)
        .def(init<af::const_ref<double, af::c_grid<2> > const&,
                  optional<double> >())
        .def(init<scitbx::sym_mat3<double> const&,
                  optional<double> >())
        .def("vectors", &w_t::vectors)
        .def("values", &w_t::values)
      ;
    }
  };

  vec3<double>
  time_eigensystem_real_symmetric(
    sym_mat3<double> const& m, std::size_t n_repetitions)
  {
    SCITBX_ASSERT(n_repetitions % 2 == 0);
    vec3<double> result(0,0,0);
    for(std::size_t i=0;i<n_repetitions/2;i++) {
      result += vec3<double>(
        eigensystem::real_symmetric<>(m).values().begin());
      result -= vec3<double>(
        eigensystem::real_symmetric<>(m).values().begin());
    }
    return result / static_cast<double>(n_repetitions);
  }

  void
  matrix_inversion_in_place_wrapper_ab(
    af::ref<double, af::c_grid<2> > const& a,
    af::ref<double, af::c_grid<2> > const& b)
  {
    if (a.accessor()[1] != a.accessor()[0]) {
      throw error("matrix_inversion_in_place: a square matrix is required.");
    }
    if (   b.accessor()[0] != 0
        && b.accessor()[1] != a.accessor()[0]) {
      throw error(
        "matrix_inversion_in_place: if a is a (n*n) matrix b must be (m*n");
    }
    if (matrix_inversion_in_place(
          a.begin(),
          static_cast<std::size_t>(a.accessor()[0]),
          b.begin(),
          static_cast<std::size_t>(b.accessor()[0])) != 0) {
      throw error("matrix is singular.");
    }
  }

  void
  matrix_inversion_in_place_wrapper_a(
    af::ref<double, af::c_grid<2> > const& a)
  {
    matrix_inversion_in_place_wrapper_ab(
      a, af::ref<double, af::c_grid<2> >(0,af::c_grid<2>(0,0)));
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    signed_phase_error_overloads, signed_phase_error, 2, 3)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    phase_error_overloads, phase_error, 2, 3)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    nearest_phase_overloads, nearest_phase, 2, 3)

  void init_module()
  {
    using namespace boost::python;

    def("floating_point_epsilon_float_get",
      &floating_point_epsilon<float>::get);
    def("floating_point_epsilon_double_get",
      &floating_point_epsilon<double>::get);

    def("erf", (double(*)(double const&)) erf);
    def("erfc", (double(*)(double const&)) erfc);
    def("erfcx", (double(*)(double const&)) erfcx);

    def("bessel_i1_over_i0", (double(*)(double const&)) bessel::i1_over_i0);
    def("bessel_inverse_i1_over_i0",
      (double(*)(double const&)) bessel::inverse_i1_over_i0);
    def("bessel_i0", (double(*)(double const&)) bessel::i0);
    def("bessel_i1", (double(*)(double const&)) bessel::i1);
    def("bessel_ln_of_i0", (double(*)(double const&)) bessel::ln_of_i0);

    eigensystem_real_symmetric_wrappers::wrap();

    wrap_gaussian();
    wrap_golay();
    wrap_minimum_covering_sphere();
    wrap_principal_axes_of_inertia();

    def("time_eigensystem_real_symmetric", time_eigensystem_real_symmetric);
    def("matrix_inversion_in_place", matrix_inversion_in_place_wrapper_ab,
      (arg_("a"), arg_("b")));
    def("matrix_inversion_in_place", matrix_inversion_in_place_wrapper_a,
      (arg_("a")));

    def("signed_phase_error",
      (double(*)(
        double const&, double const&, bool))
          math::signed_phase_error,
      signed_phase_error_overloads(
        (arg_("phi1"), arg_("phi2"), arg_("deg")=false)));
    def("signed_phase_error",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::signed_phase_error,
      signed_phase_error_overloads(
        (arg_("phi1"), arg_("phi2"), arg_("deg")=false)));
    def("phase_error",
      (double(*)(
        double const&, double const&, bool))
          math::phase_error,
      phase_error_overloads(
        (arg_("phi1"), arg_("phi2"), arg_("deg")=false)));
    def("phase_error",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::phase_error,
      phase_error_overloads(
        (arg_("phi1"), arg_("phi2"), arg_("deg")=false)));
    def("nearest_phase",
      (double(*)(
        double const&, double const&, bool))
          math::nearest_phase,
      nearest_phase_overloads(
        (arg_("reference"), arg_("other"), arg_("deg")=false)));
    def("nearest_phase",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::nearest_phase,
      nearest_phase_overloads(
        (arg_("reference"), arg_("other"), arg_("deg")=false)));
  }

}}}} // namespace scitbx::math::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_math_ext)
{
  scitbx::math::boost_python::init_module();
}
