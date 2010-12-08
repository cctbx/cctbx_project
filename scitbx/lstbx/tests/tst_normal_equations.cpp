#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/array_family/initialiser.h>
#include <scitbx/array_family/simple_io.h>
#include <iostream>
#include <iomanip>

using namespace scitbx;
using scitbx::fn::approx_equal;

template <template<typename> class NormalEquationsType>
struct linear_polynomial_fit
{
  double noise;
  static int const n_data = 10, n_params=3;
  typedef NormalEquationsType<double> normal_eqns_t;
  typedef typename normal_eqns_t::symmetric_matrix_owning_ref_t
          symmetric_matrix_owning_ref_t;
  typedef typename normal_eqns_t::upper_diagonal_matrix_owning_ref_t
          upper_diagonal_matrix_owning_ref_t;
  typedef typename normal_eqns_t::vector_owning_ref_t
          vector_owning_ref_t;
  normal_eqns_t ls;

  linear_polynomial_fit(double noise_)
    : noise(noise_), ls(n_params, false)
  {}

  void compute(double a, double b, double c) {
    double sign = +1;
    double w = 1.;
    for (int i=0; i<n_data; ++i) {
      double t = i * 1./n_data;
      double yo = 2*t*t*(1-t) + noise*sign;
      sign = -sign;
      double t2=t*t, t3=t*t2;
      double yc = -t3 + a*t2 + b*t + c;
      double grad_yc[n_params] = { t2, t, 1 };
      ls.add_equation(yc, grad_yc, yo, w);
    }
    ls.finalise();
  }

  static void exercise()
  {
    // the optimal scale factor is 2
    double const tol=1e-15;
    bool exception_thrown;
    {
      linear_polynomial_fit fit(1e-5);
      fit.compute(0.5, 0.3, 0.2);
      SCITBX_ASSERT(
        approx_equal(fit.ls.objective(), 0.039642707534326034, 5e-14));
      SCITBX_ASSERT(
        approx_equal(fit.ls.optimal_scale_factor(),
                     0.6148971786833856, 5e-14));
      lstbx::normal_equations::linear_ls<double>
      &normal_eqns = fit.ls.step_equations();
      SCITBX_ASSERT(!normal_eqns.solved());
      try {
        normal_eqns.cholesky_factor();
        exception_thrown = false;
      }
      catch (scitbx::error) {
        exception_thrown = true;
      }
      SCITBX_ASSERT(exception_thrown);
      try {
        normal_eqns.solution();
        exception_thrown = false;
      }
      catch (scitbx::error) {
        exception_thrown = true;
      }
      SCITBX_ASSERT(exception_thrown);
      // C.f. Mathematica notebook tst_normal_equations for how
      // these reference values have been obtained.
      symmetric_matrix_owning_ref_t a = normal_eqns.normal_matrix();
      vector_owning_ref_t b = normal_eqns.right_hand_side();
      symmetric_matrix_owning_ref_t a0(3);
      af::init(a0) =
        0.371944193675858, 0.39066546997866547, 0.10797294655500618,
                           0.41859250354804045, 0.08077629438075473,
                                                0.19767268057900367;
      vector_owning_ref_t b0(3);
      af::init(b0) =
        0.12149917297914861, 0.13803759252793774, -0.025190641142579157;
      SCITBX_ASSERT(a.all_approx_equal(a0, 5e-14));
      SCITBX_ASSERT(b.all_approx_equal(b0, 5e-14));

      normal_eqns.solve();
      SCITBX_ASSERT(normal_eqns.solved());
      try {
        normal_eqns.normal_matrix();
        exception_thrown = false;
      }
      catch (scitbx::error) {
        exception_thrown = true;
      }
      SCITBX_ASSERT(exception_thrown);
      try {
        normal_eqns.right_hand_side();
        exception_thrown = false;
      }
      catch (scitbx::error) {
        exception_thrown = true;
      }
      SCITBX_ASSERT(exception_thrown);
      upper_diagonal_matrix_owning_ref_t u0(3), u=normal_eqns.cholesky_factor();
      af::init(u0) =
        0.6098722765266986, 0.6405693208478925 ,  0.1770418999366983 ,
                            0.09090351333425013, -0.3589664912436558 ,
                                                  0.19357661121640218;
      SCITBX_ASSERT(u.all_approx_equal(u0, 5e-14));
      vector_owning_ref_t s0(3), s=normal_eqns.solution();
      af::init(s0) =
        1.2878697604109028, -0.7727798877778043, -0.5151113342942297;
      SCITBX_ASSERT(s.all_approx_equal(s0, 1e-12));
    }
  }
};

struct linear_combination_fit
{
  static int const n=4, p=3, n_data=10;
  double noise;
  typedef lstbx::normal_equations::separable_non_linear_ls<double>
          normal_eqns_t;
  typedef normal_eqns_t::symmetric_matrix_owning_ref_t
          symmetric_matrix_owning_ref_t;
  typedef normal_eqns_t::vector_owning_ref_t
          vector_owning_ref_t;
  typedef normal_eqns_t::vector_t vector_t;
  normal_eqns_t ls;

  linear_combination_fit(double noise_)
    : noise(noise_), ls(n, p)
  {}

  void compute(double x1, double x2, double x3, double x4) {
    double sign = +1;
    double w = 1.;
    for (int i=0; i<n_data; ++i) {
      double t0 = i * 1./n_data;
      af::tiny<double, p> t(t0, 1-t0, -t0) ;
      double yo = noise*sign;
      for (int k=0; k<p; ++k) yo += t[k]*(1 + t[k]*(-1 + 2*t[k]));
      sign = -sign;

      af::tiny<double, p> yc;;
      af::versa<double, af::mat_grid> grad_yc_(af::mat_grid(n,p));
      af::ref<double, af::mat_grid> grad_yc = grad_yc_.ref();
      for (int k=0; k<p; ++k) {
        double t1 = t[k];
        double t2 = t1*t1;
        double t3 = t1*t2;
        double t4 = t2*t2;
        yc[k] = x4*t4 + x3*t3 + x2*t2 + x1*t1;
        grad_yc(0, k) += t1;
        grad_yc(1, k) += t2;
        grad_yc(2, k) += t3;
        grad_yc(3, k) += t4;
      }
      ls.add_equation(yc.ref(), grad_yc, yo, w);
    }
  }

  static void exercise() {
    linear_combination_fit fit(1e-5);
    fit.compute(1, 2, 3, 4);
    vector_t a_star = fit.ls.optimal_linear_coefficients();
    vector_t a_star_0(p);
    af::init(a_star_0) = -0.210323669833525,
                          0.208033764252907,
                         -0.0543166780808711;
    SCITBX_ASSERT(a_star.all_approx_equal(a_star_0, 1e-14));
    lstbx::normal_equations::linear_ls<double> normal_eqns = fit.ls.equations();
    symmetric_matrix_owning_ref_t a = normal_eqns.normal_matrix();
    vector_owning_ref_t b = normal_eqns.right_hand_side();

    vector_owning_ref_t b0(n);
    af::init(b0) = 0.00542544409277865,   0.00539923674666925,
                   0.000326266092452357, -0.00430067896586856;
    SCITBX_ASSERT(b.all_approx_equal(b0, 1e-14));

    symmetric_matrix_owning_ref_t a0(n);
    af::init(a0) =
      0.00300136362055691, 0.00227988574887292, 0.000520356169835869,
                                                          -0.000924189883757925,

                           0.00222887407045411, 0.000625283216528185,
                                                          -0.000803561698174109,

                                                0.0000743095490471244,
                                                          -0.000416896289395314,

                                                           -0.000129669204394118;
    SCITBX_ASSERT(a.all_approx_equal(a0, 1e-14))(a.const_ref());
  }
};

int main() {
  linear_polynomial_fit<
    lstbx::normal_equations
    ::non_linear_ls_with_separable_scale_factor>::exercise();
  linear_combination_fit::exercise();
  std::cout << "OK\n";
  return 0;
}
