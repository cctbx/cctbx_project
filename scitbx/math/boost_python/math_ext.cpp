#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/math/erf.h>
#include <scitbx/math/bessel.h>
#include <scitbx/math/gamma.h>
#include <scitbx/math/chebyshev.h>
#include <scitbx/math/lambertw.h>
#include <scitbx/math/superpose.h>
#include <scitbx/math/phase_error.h>
#include <scitbx/math/resample.h>
#include <scitbx/math/halton.h>
#include <scitbx/math/utils.h>
#include <scitbx/math/euler_angles.h>
#include <scitbx/math/gcd.h>
#include <boost/rational.hpp> // for boost::gcd

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

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
  void wrap_r3_rotation();
  void wrap_resample();
  void wrap_quadrature();
  void wrap_unimodular_generator();
  void wrap_halton();
  void wrap_least_squares_plane();
  void wrap_continued_fraction();
  void wrap_numeric_limits();

namespace {

  int
  time_gcd_int_boost_math(
    int n)
  {
    int result = 0;
    for(int a=0;a<n;a++) {
      for(int b=0;b<n;b++) {
        int c = boost::gcd(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }

  int
  time_gcd_int_simple(
    int n)
  {
    int result = 0;
    for(int a=0;a<n;a++) {
      for(int b=0;b<n;b++) {
        int c = gcd_int_simple(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }

#if defined(SCITBX_MATH_GCD_USING_ASM)
  int
  time_gcd_int_asm(
    int n)
  {
    int result = 0;
    for(int a=0;a<n;a++) {
      for(int b=0;b<n;b++) {
        int c = gcd_int_asm(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }
#endif

  mat3<double>
  superpose_kearsley_rotation(
    af::const_ref<vec3<double> > const& reference_sites,
    af::const_ref<vec3<double> > const& other_sites)
  {
    return superpose::superposition<>::kearsley_rotation(
      reference_sites, other_sites);
  }

  mat3< double >
  euler_angles_xyz_matrix(
      const double& ax,
      const double& ay,
      const double& az )
  {
    return euler_angles::xyz_matrix( ax, ay, az );
  }

  vec3< double >
  euler_angles_xyz_angles(
      const mat3< double >& m,
      const double& eps = 1e-12 )
  {
    return euler_angles::xyz_angles( m, eps );
  }

  mat3< double >
  euler_angles_yzx_matrix(
      const double& ay,
      const double& az,
      const double& ax )
  {
    return euler_angles::yzx_matrix( ay, az, ax );
  }

  vec3< double >
  euler_angles_yzx_angles(
      const mat3< double >& m,
      const double& eps = 1e-12 )
  {
    return euler_angles::yzx_angles( m, eps );
  }

  mat3< double >
  euler_angles_zyz_matrix(
      const double& az1,
      const double& ay,
      const double& az3 )
  {
    return euler_angles::zyz_matrix( az1, ay, az3 );
  }

  vec3< double >
  euler_angles_zyz_angles(
      const mat3< double >& m,
      const double& eps = 1e-12 )
  {
    return euler_angles::zyz_angles( m, eps );
  }

  void init_module()
  {
    using namespace boost::python;

    def("time_gcd_int_boost_math", time_gcd_int_boost_math);
    def("gcd_int_simple", gcd_int_simple, (arg("a"), arg("b")));
    def("time_gcd_int_simple", time_gcd_int_simple);
#if defined(SCITBX_MATH_GCD_USING_ASM)
    def("gcd_int_asm", gcd_int_asm, (arg("a"), arg("b")));
    def("time_gcd_int_asm", time_gcd_int_asm);
#endif

    def("floating_point_epsilon_float_get",
      &floating_point_epsilon<float>::get);
    def("floating_point_epsilon_double_get",
      &floating_point_epsilon<double>::get);

    def("erf", (double(*)(double const&)) erf);
    def("erf", (scitbx::af::shared<double>(*)(scitbx::af::const_ref<double> const&)) erf);

    def("erfc", (double(*)(double const&)) erfc);
    def("erfcx", (double(*)(double const&)) erfcx);

    def("bessel_i1_over_i0", (double(*)(double const&)) bessel::i1_over_i0);
    def("bessel_i1_over_i0", (scitbx::af::shared<double>(*)(scitbx::af::const_ref<double> const&)) bessel::i1_over_i0);
    def("bessel_inverse_i1_over_i0",
      (double(*)(double const&)) bessel::inverse_i1_over_i0);
    def("inverse_bessel_i1_over_i0", (scitbx::af::shared<double>(*)(
         scitbx::af::const_ref<double> const&)) bessel::inverse_i1_over_i0);
    def("bessel_i0", (double(*)(double const&)) bessel::i0);
    def("bessel_i1", (double(*)(double const&)) bessel::i1);
    def("bessel_ln_of_i0", (double(*)(double const&)) bessel::ln_of_i0);
    def("ei1", (double(*)(double const&)) bessel::ei1);
    def("ei0", (double(*)(double const&)) bessel::ei0);



    def("gamma_complete", (double(*)(double const&, bool))
      gamma::complete, (
        arg_("x"),
        arg_("minimax")=true));
    def("gamma_incomplete", (double(*)(double const&,
                                       double const&,
                                       unsigned))
      gamma::incomplete, (arg_("a"),
        arg_("x"),
        arg_("max_iterations")=500));
    def("gamma_incomplete_complement",(double(*)(double const&,
                                                 double const&,
                                                 unsigned))
      gamma::incomplete_complement, (
        arg_("a"),
        arg_("x"),
        arg_("max_iterations")=500));
    def("exponential_integral_e1z", (double(*)(double const&))
         gamma::exponential_integral_e1z );



    def("lambertw", (double(*)(double const&, unsigned))
      lambertw, (
        arg_("x"),
        arg_("max_iterations")=100));

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
    wrap_r3_rotation();
    wrap_resample();
    wrap_quadrature();
    wrap_unimodular_generator();
    wrap_halton();
    wrap_least_squares_plane();
    wrap_continued_fraction();
    wrap_numeric_limits();

    def("superpose_kearsley_rotation", superpose_kearsley_rotation, (
      arg_("reference_sites"), arg_("other_sites")));

    def( "euler_angles_xyz_matrix", euler_angles_xyz_matrix, (
      arg_("ax"), arg_("ay"), arg_("az")));
    def( "euler_angles_xyz_angles", euler_angles_xyz_angles, (
      arg_("m"), arg_("eps")=1e-12));

    def( "euler_angles_yzx_matrix", euler_angles_yzx_matrix, (
      arg_("ay"), arg_("az"), arg_("ax")));
    def( "euler_angles_yzx_angles", euler_angles_yzx_angles, (
      arg_("m"), arg_("eps")=1e-12));

    def( "euler_angles_zyz_matrix", euler_angles_zyz_matrix, (
      arg_("az1"), arg_("ay"), arg_("az3")));
    def( "euler_angles_zyz_angles", euler_angles_zyz_angles, (
      arg_("m"), arg_("eps")=1e-12));

    def("signed_phase_error",
      (double(*)(
        double const&, double const&, bool))
          math::signed_phase_error, (
            arg_("phi1"), arg_("phi2"), arg_("deg")=false));
    def("signed_phase_error",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::signed_phase_error, (
            arg_("phi1"), arg_("phi2"), arg_("deg")=false));
    def("phase_error",
      (double(*)(
        double const&, double const&, bool))
          math::phase_error, (
            arg_("phi1"), arg_("phi2"), arg_("deg")=false));
    def("phase_error",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::phase_error, (
            arg_("phi1"), arg_("phi2"), arg_("deg")=false));
    def("nearest_phase",
      (double(*)(
        double const&, double const&, bool))
          math::nearest_phase, (
            arg_("reference"), arg_("other"), arg_("deg")=false));
    def("nearest_phase",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::nearest_phase, (
            arg_("reference"), arg_("other"), arg_("deg")=false));
    def("divmod", math::divmod);
  }

}}}} // namespace scitbx::math::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_math_ext)
{
  scitbx::math::boost_python::init_module();
}
