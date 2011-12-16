#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/bessel.h>
#include <scitbx/math/chebyshev.h>
#include <scitbx/math/dihedral.h>
#include <scitbx/math/erf.h>
#include <scitbx/math/euler_angles.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/math/gamma.h>
#include <scitbx/math/gcd.h>
#include <scitbx/math/halton.h>
#include <scitbx/math/lambertw.h>
#include <scitbx/math/phase_error.h>
#include <scitbx/math/resample.h>
#include <scitbx/math/superpose.h>
#include <scitbx/math/utils.h>
#include <boost/rational.hpp> // for boost::gcd
#include <scitbx/math/approx_equal.h>
#include <scitbx/math/orthonormal_basis.h>
#include <scitbx/math/gaussian_fit_1d_analytical.h>
#include <scitbx/math/cubic_equation.h>
#include <scitbx/math/distance_difference.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

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
  void wrap_distributions();
  void wrap_exp_functions();
  void wrap_zernike();
  void wrap_zernike_mom();
  void wrap_2d_zernike_mom();
  void wrap_weighted_covariance();
  void wrap_dmatrix();
  void wrap_correlation();
  void wrap_interpolation();

namespace {

  int
  time_gcd_int_boost(
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

  long
  time_gcd_long_boost(
    long n)
  {
    long result = 0;
    for(long a=0;a<n;a++) {
      for(long b=0;b<n;b++) {
        long c = boost::gcd(a, b);
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

  long
  time_gcd_long_simple(
    long n)
  {
    long result = 0;
    for(long a=0;a<n;a++) {
      for(long b=0;b<n;b++) {
        long c = gcd_long_simple(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }

  long
  time_gcd_unsigned_long_binary(
    unsigned long n)
  {
    unsigned long result = 0;
    for(unsigned long a=0;a<n;a++) {
      for(unsigned long b=0;b<n;b++) {
        unsigned long c = gcd_unsigned_long_binary(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }

  long
  time_gcd_long_binary(
    long n)
  {
    long result = 0;
    for(long a=0;a<n;a++) {
      for(long b=0;b<n;b++) {
        long c = gcd_long_binary(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }

#if defined(SCITBX_MATH_GCD_USING_ASM)
  int
  time_gcd_int32_asm(
    int n)
  {
    int result = 0;
    for(int a=0;a<n;a++) {
      for(int b=0;b<n;b++) {
        int c = gcd_int32_asm(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }

# if defined(__x86_64__)
  long
  time_gcd_int64_asm(
    long n)
  {
    long result = 0;
    for(long a=0;a<n;a++) {
      for(long b=0;b<n;b++) {
        long c = gcd_int64_asm(a, b);
        if (result < c) result = c;
      }
    }
    return result;
  }
# endif
#endif

  template <typename SitesType>
  boost::optional<double>
  dihedral_angle(
    SitesType const& sites,
    bool deg)
  {
    return dihedral(sites).angle(deg);
  }

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

  template <typename T>
  struct approx_equal_relatively_wrapper
  {
    typedef math::approx_equal_relatively<T> wt;
    typedef T arg_t;
    typedef typename wt::amplitude_type amplitude_t;
    static bool form_1(arg_t x, arg_t y, amplitude_t relative_error) {
      wt p(relative_error);
      return p(x, y);
    }
    static bool form_2(arg_t x, arg_t y, amplitude_t relative_error,
                       amplitude_t near_zero_threshold)
    {
      wt p(relative_error, near_zero_threshold);
      return p(x, y);
    }

    static void wrap() {
      using namespace boost::python;
      def("approx_equal_relatively", form_1,
          (arg("x"), arg("y"), arg("relative_error")));
      def("approx_equal_relatively", form_2,
          (arg("x"), arg("y"), arg("relative_error"),
           arg("near_zero_threshold")));
    }
  };

  void init_module()
  {
    using namespace boost::python;

    def("time_gcd_int_boost", time_gcd_int_boost);
    def("time_gcd_long_boost", time_gcd_long_boost);
    def("gcd_int_simple", gcd_int_simple, (arg("a"), arg("b")));
    def("time_gcd_int_simple", time_gcd_int_simple);
    def("gcd_long_simple", gcd_long_simple, (arg("a"), arg("b")));
    def("time_gcd_long_simple", time_gcd_long_simple);
    def("time_gcd_unsigned_long_binary", time_gcd_unsigned_long_binary);
    def("gcd_long_binary", gcd_long_binary, (arg("a"), arg("b")));
    def("time_gcd_long_binary", time_gcd_long_binary);
#if defined(SCITBX_MATH_GCD_USING_ASM)
    def("gcd_int32_asm", gcd_int32_asm, (arg("a"), arg("b")));
    def("time_gcd_int32_asm", time_gcd_int32_asm);
# if defined(__x86_64__)
    def("gcd_int64_asm", gcd_int64_asm, (arg("a"), arg("b")));
    def("time_gcd_int64_asm", time_gcd_int64_asm);
# endif
#endif

    def("floating_point_epsilon_float_get",
      &floating_point_epsilon<float>::get);
    def("floating_point_epsilon_double_get",
      &floating_point_epsilon<double>::get);

    def("erf", (double(*)(double const&)) erf);
    def("erf",
      (scitbx::af::shared<double>(*)(
        scitbx::af::const_ref<double> const&)) erf);
    def("erfc", (double(*)(double const&)) erfc);
    def("erfcx", (double(*)(double const&)) erfcx);

    def("bessel_i1_over_i0", (double(*)(double const&)) bessel::i1_over_i0);
    def("bessel_i1_over_i0",
      (scitbx::af::shared<double>(*)(scitbx::af::const_ref<double> const&))
        bessel::i1_over_i0);
    def("bessel_inverse_i1_over_i0",
      (double(*)(double const&)) bessel::inverse_i1_over_i0);
    def("inverse_bessel_i1_over_i0", (scitbx::af::shared<double>(*)(
         scitbx::af::const_ref<double> const&)) bessel::inverse_i1_over_i0);
    def("bessel_i0", (double(*)(double const&)) bessel::i0);
    def("bessel_i1", (double(*)(double const&)) bessel::i1);
    def("bessel_ln_of_i0", (double(*)(double const&)) bessel::ln_of_i0);
    def("ei1", (double(*)(double const&)) bessel::ei1);
    def("ei0", (double(*)(double const&)) bessel::ei0);

    typedef return_value_policy<return_by_value> rbv;
    namespace smg=scitbx::math::gaussian_fit_1d_analytical;
    class_<smg::compute<> >("gaussian_fit_1d_analytical")
      .def(init<
           af::const_ref<double> const&,
           af::const_ref<double> const& >(
             (arg("x"),
              arg("y"))))
      .add_property("a", make_getter(&smg::compute<>::a, rbv()))
      .add_property("b", make_getter(&smg::compute<>::b, rbv()))
    ;

    //typedef return_value_policy<return_by_value> rbv;
    namespace cueq=scitbx::math::cubic_equation;
    class_<cueq::real<> >("cubic_equation_real")
      .def(init<
           double const&,
           double const&,
           double const&,
           double const& >(
             (arg("a"),
              arg("b"),
              arg("c"),
              arg("d"))))
      .def("residual", &cueq::real<>::residual)
      .add_property("x", make_getter(&cueq::real<>::x, rbv()))
      .add_property("A", make_getter(&cueq::real<>::A,  rbv()))
      .add_property("B", make_getter(&cueq::real<>::B,  rbv()))
      .add_property("D", make_getter(&cueq::real<>::D,  rbv()))
    ;

#if defined(SCITBX_MATH_BESSEL_HAS_SPHERICAL)
    def("spherical_bessel",
      (double(*)(int const&, double const&))
        bessel::spherical_bessel);
    def("spherical_bessel_array",
      (scitbx::af::shared< double> (*)(
        int const&, scitbx::af::shared<double> const&))
          bessel::spherical_bessel_array);
    def("bessel_J",
      (double(*)(int const&, double const&))
        bessel::bessel_J);
    def("bessel_J_array",
      (scitbx::af::shared< double> (*)(
        int const&, scitbx::af::shared<double> const&))
          bessel::bessel_J_array);

#endif

    def("gamma_complete", (double(*)(double const&, bool))
      gamma::complete, (
        arg("x"),
        arg("minimax")=true));
    def("gamma_incomplete", (double(*)(double const&,
                                       double const&,
                                       unsigned))
      gamma::incomplete, (arg("a"),
        arg("x"),
        arg("max_iterations")=500));
    def("gamma_incomplete_complement",(double(*)(double const&,
                                                 double const&,
                                                 unsigned))
      gamma::incomplete_complement, (
        arg("a"),
        arg("x"),
        arg("max_iterations")=500));
    def("exponential_integral_e1z", (double(*)(double const&))
         gamma::exponential_integral_e1z );

    def("lambertw", (double(*)(double const&, unsigned))
      lambertw, (
        arg("x"),
        arg("max_iterations")=100));

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
    wrap_distributions();
    wrap_exp_functions();
    wrap_zernike();
    wrap_zernike_mom();
    wrap_2d_zernike_mom();
    wrap_weighted_covariance();
    wrap_dmatrix();
    wrap_correlation();
    wrap_interpolation();

    def("superpose_kearsley_rotation", superpose_kearsley_rotation, (
      arg("reference_sites"), arg("other_sites")));

    def("dihedral_angle",
      dihedral_angle<af::tiny<vec3<double>, 4> >, (
        arg("sites"), arg("deg")=false));
    def("dihedral_angle",
      dihedral_angle<af::const_ref<vec3<double> > >, (
        arg("sites"), arg("deg")=false));

    def( "euler_angles_xyz_matrix", euler_angles_xyz_matrix, (
      arg("ax"), arg("ay"), arg("az")));
    def( "euler_angles_xyz_angles", euler_angles_xyz_angles, (
      arg("m"), arg("eps")=1e-12));

    def( "euler_angles_yzx_matrix", euler_angles_yzx_matrix, (
      arg("ay"), arg("az"), arg("ax")));
    def( "euler_angles_yzx_angles", euler_angles_yzx_angles, (
      arg("m"), arg("eps")=1e-12));

    def( "euler_angles_zyz_matrix", euler_angles_zyz_matrix, (
      arg("az1"), arg("ay"), arg("az3")));
    def( "euler_angles_zyz_angles", euler_angles_zyz_angles, (
      arg("m"), arg("eps")=1e-12));

    def("signed_phase_error",
      (double(*)(
        double const&, double const&, bool))
          math::signed_phase_error, (
            arg("phi1"), arg("phi2"), arg("deg")=false));
    def("signed_phase_error",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::signed_phase_error, (
            arg("phi1"), arg("phi2"), arg("deg")=false));
    def("phase_error",
      (double(*)(
        double const&, double const&, bool))
          math::phase_error, (
            arg("phi1"), arg("phi2"), arg("deg")=false));
    def("phase_error",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::phase_error, (
            arg("phi1"), arg("phi2"), arg("deg")=false));
    def("nearest_phase",
      (double(*)(
        double const&, double const&, bool))
          math::nearest_phase, (
            arg("reference"), arg("other"), arg("deg")=false));
    def("nearest_phase",
      (af::shared<double>(*)(
        af::const_ref<double> const&, af::const_ref<double> const&, bool))
          math::nearest_phase, (
            arg("reference"), arg("other"), arg("deg")=false));
    def("divmod", math::divmod);
    approx_equal_relatively_wrapper<double>::wrap();
    approx_equal_relatively_wrapper<std::complex<double> >::wrap();
    {
      af::tiny<vec3<double>, 3> (*f1)(vec3<double> const &, vec3<double> const &,
                                      bool) = &orthonormal_basis;
      af::tiny<vec3<double>, 3> (*f2)(vec3<double> const &, int,
                                      vec3<double> const &, int,
                                      bool) = &orthonormal_basis;
      def("orthonormal_basis", f1, (arg("v0"), arg("v1"),
                                    arg("right_handed")=true));
      def("orthonormal_basis", f2, (arg("v0"), arg("axis_index_1"),
                                    arg("v1"), arg("axis_index_2"),
                                    arg("right_handed")=true));
    }
    def("distance_difference_matrix", distance_difference_matrix<double>, (
      arg("sites1"), arg("sites2")));
  }

}}}} // namespace scitbx::math::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_math_ext)
{
  scitbx::math::boost_python::init_module();
}
