#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/math/erf.h>
#include <scitbx/math/bessel.h>
#include <scitbx/math/gamma.h>
#include <scitbx/math/chebyshev.h>
#include <scitbx/math/lambertw.h>
#include <scitbx/math/eigensystem.h>
#include <scitbx/math/phase_error.h>
#include <scitbx/math/resample.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace math {
namespace boost_python {

  void wrap_basic_statistics();
  void wrap_gaussian();
  void wrap_golay();
  void wrap_minimum_covering_sphere();
  void wrap_principal_axes_of_inertia();
  void wrap_row_echelon();
  void wrap_tensor_rank_2();
  void wrap_icosahedron();
  void wrap_chebyshev_base();
  void wrap_chebyshev_polynome();
  void wrap_chebyshev_fitter();
  void wrap_chebyshev_lsq();
  void wrap_slatec();
  void wrap_line_search();
  void wrap_non_parametric_bootstrap();
  void wrap_smooth_bootstrap();

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

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    gamma_complete_overloads,
    gamma::complete, 1, 2)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    gamma_incomplete_overloads,
    gamma::incomplete, 2, 3)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    gamma_incomplete_complement_overloads,
    gamma::incomplete_complement, 2, 3)

  BOOST_PYTHON_FUNCTION_OVERLOADS(lambertw_overloads, lambertw, 1, 2)

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
    def("erf", (scitbx::af::shared<double>(*)(scitbx::af::const_ref<double> const&)) erf);

    def("erfc", (double(*)(double const&)) erfc);
    def("erfcx", (double(*)(double const&)) erfcx);

    def("bessel_i1_over_i0", (double(*)(double const&)) bessel::i1_over_i0);
    def("bessel_inverse_i1_over_i0",
      (double(*)(double const&)) bessel::inverse_i1_over_i0);
    def("bessel_i0", (double(*)(double const&)) bessel::i0);
    def("bessel_i1", (double(*)(double const&)) bessel::i1);
    def("bessel_ln_of_i0", (double(*)(double const&)) bessel::ln_of_i0);

    def("gamma_complete", (double(*)(double const&, bool))
        gamma::complete,
        gamma_complete_overloads( (arg_("x"),
                                   arg_("minimax")=true)));
    def("gamma_incomplete", (double(*)(double const&,
                                       double const&,
                                       unsigned))
        gamma::incomplete,
        gamma_incomplete_overloads( (arg_("a"),
                                     arg_("x"),
                                     arg_("max_iterations")=500 )));
    def("gamma_incomplete_complement",(double(*)(double const&,
                                                 double const&,
                                                 unsigned))
        gamma::incomplete_complement,
        gamma_incomplete_complement_overloads( (arg_("a"),
                                                arg_("x"),
                                                arg_("max_iterations")=500 )));




    def("lambertw", (double(*)(double const&, unsigned)) lambertw,
      lambertw_overloads(
        (arg_("x"), arg_("max_iterations")=100)));

    eigensystem_real_symmetric_wrappers::wrap();

    wrap_basic_statistics();
    wrap_gaussian();
    wrap_golay();
    wrap_minimum_covering_sphere();
    wrap_principal_axes_of_inertia();
    wrap_row_echelon();
    wrap_tensor_rank_2();
    wrap_icosahedron();
    wrap_chebyshev_base();
    wrap_chebyshev_polynome();
    wrap_chebyshev_fitter();
    wrap_chebyshev_lsq();
    wrap_slatec();
    wrap_line_search();
    // resampling
    wrap_non_parametric_bootstrap();
    wrap_smooth_bootstrap();

    def("time_eigensystem_real_symmetric", time_eigensystem_real_symmetric);

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
