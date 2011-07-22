#include <cmath>
#include <scitbx/sym_mat3.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/array_family/simple_io.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int /*argc*/, char* /*argv*/[])
{
  {
    sym_mat3<int> va;
    check_true(__LINE__, va.size() == 6);
    sym_mat3<int> vb(0,0,0,0,0,0);
    for(int j=0;j<6;j++) check_true(__LINE__, vb[j] == 0);
    sym_mat3<int> vc(af::tiny_plain<int,6>(0,0,0,0,0,0));
    for(int j=0;j<6;j++) check_true(__LINE__, vc[j] == 0);
    int id[] = {1,1,1,1,1,1};
    sym_mat3<int> vd(id);
    for(int j=0;j<6;j++) check_true(__LINE__, vd[j] == 1);
    mat3<int> ve(sym_mat3<int>(1,2,3,4,5,6));
    check_true(__LINE__, ve == mat3<int>(1,4,5, 4,2,6, 5,6,3));
  }
  {
    sym_mat3<int> a; a.fill(0);
    sym_mat3<int> b; b.fill(0);
    sym_mat3<int> c; c.fill(1);
    check_true(__LINE__, a == b);
    check_false(__LINE__, a == c);
    check_true(__LINE__, a == 0);
    check_false(__LINE__, a == 1);
    check_true(__LINE__, 0 == a);
    check_false(__LINE__, 1 == a);
    check_false(__LINE__, a != b);
    check_true(__LINE__, a != c);
    check_false(__LINE__, a != 0);
    check_true(__LINE__, a != 1);
    check_false(__LINE__, 0 != a);
    check_true(__LINE__, 1 != a);
    check_true(__LINE__, sym_mat3<int>(1) == sym_mat3<int>(1,1,1,0,0,0));
    check_true(__LINE__, sym_mat3<int>(vec3<int>(1,2,3)) == sym_mat3<int>(
      1,2,3,0,0,0));
  }
  {
    sym_mat3<int> a(1,2,3,4,5,6);
    sym_mat3<int> b(4,5,6,7,8,9);
    sym_mat3<int> c(5,7,9,11,13,15);
    sym_mat3<int> d(2,4,6,8,10,12);
    check_true(__LINE__, a(0,0) == 1);
    check_true(__LINE__, a(0,1) == 4);
    check_true(__LINE__, a(0,2) == 5);
    check_true(__LINE__, a(1,0) == 4);
    check_true(__LINE__, a(1,1) == 2);
    check_true(__LINE__, a(1,2) == 6);
    check_true(__LINE__, a(2,0) == 5);
    check_true(__LINE__, a(2,1) == 6);
    check_true(__LINE__, a(2,2) == 3);
    check_true(__LINE__, a + b == c);
    check_true(__LINE__, a + 3 == b);
    check_true(__LINE__, 3 + a == b);
    check_true(__LINE__, a - b == -3);
    check_true(__LINE__, b - 3 == a);
    check_true(__LINE__, c - b == a);
    check_true(__LINE__, a * b == mat3<int>(
      72,72,74,78,92,86,86,92,112));
    check_true(__LINE__, b * a == mat3<int>(
      72,78,86,72,92,92,74,86,112));
    {
      mat3<int> f(1,2,3, 3,-4,5, 4,5,6);
      check_true(__LINE__, f * a == f * mat3<int>(a));
      check_true(__LINE__, a * f == mat3<int>(a) * f);
    }
    check_true(__LINE__, a * vec3<int>(1,2,3) == vec3<int>(
      24,26,26));
    check_true(__LINE__, std::fabs(af::max(
      (b * vec3<double>(3,5,7) - vec3<double>(103,109,111)).ref())) < 1.e-6);
    check_true(__LINE__, vec3<int>(1,2,3) * a == vec3<int>(
      24,26,26));
    check_true(__LINE__, std::fabs(af::max(
      (vec3<double>(3,5,7) * b - vec3<double>(103,109,111)).ref())) < 1.e-6);
    check_true(__LINE__, a * 2 == d);
    check_true(__LINE__, 2 * a == d);
    check_true(__LINE__, d / 2 == a);
    check_true(__LINE__, 5040 / d == sym_mat3<int>(
      2520,1260,840,630,504,420));
    check_true(__LINE__, a % 2 == sym_mat3<int>(1,0,1,0,1,0));
    check_true(__LINE__, 4 % d == sym_mat3<int>(0,0,4,4,4,4));
    {
      sym_mat3<int> t(a);
      check_true(__LINE__, (t += b) == c);
    }
    {
      sym_mat3<int> t(a);
      check_true(__LINE__, (t += 3) == b);
    }
    {
      sym_mat3<int> t(a);
      check_true(__LINE__, (t -= b) == -3);
    }
    {
      sym_mat3<int> t(b);
      check_true(__LINE__, (t -= 3) == a);
    }
    {
      sym_mat3<int> t(a);
      check_true(__LINE__, (t *= 2) == d);
    }
    {
      sym_mat3<int> t(d);
      check_true(__LINE__, (t /= 2) == a);
    }
    {
      sym_mat3<int> t(a);
      check_true(__LINE__, (t %= 2) == sym_mat3<int>(1,0,1,0,1,0));
    }
    check_true(__LINE__, -a == sym_mat3<int>(-1,-2,-3,-4,-5,-6));
    check_true(__LINE__, +a == a);
    check_true(__LINE__, a.diagonal() == vec3<int>(1,2,3));
    check_true(__LINE__, a.transpose() == a);
    check_true(__LINE__, a.trace() == 6);
    check_true(__LINE__, a.determinant() == 112);
    check_true(__LINE__, sym_mat3<int>(1).determinant() == 1);
    check_true(__LINE__, sym_mat3<int>(2).determinant() == 8);
    check_true(__LINE__, a.co_factor_matrix_transposed() == sym_mat3<int>(
      -30, -22, -14, 18, 14, 14));
    sym_mat3<double> fa(a);
    sym_mat3<double> fai = fa.inverse();
    check_true(__LINE__, std::fabs(af::max(
      ((fai * fa) - mat3<double>(1)).ref())) < 1.e-6);
    check_true(__LINE__, a.dot(b) == 276);
    check_true(__LINE__, b.dot(a) == 276);
  }
  {
    mat3<double> t;
    t.set_row(0, vec3<int>(1,2,3));
    t.set_row(1, vec3<int>(2,5,6));
    t.set_row(2, vec3<int>(3,6,9));
    check_true(__LINE__, sym_mat3<double>(t).determinant() == t.determinant());
    t[1] = 2.01;
    sym_mat3<double> s(t, 0.1);
    t[1] = 2.005;
    t[3] = 2.005;
    check_true(__LINE__, std::fabs(s.determinant() - t.determinant()) < 1.e-6);
  }
  {
    using scitbx::matrix::matrix_vector;
    sym_mat3<int> a(3,4,9, 2,1,7);
    mat3<int> c(1,2,3, 3,-4,5, 4,5,6);
    sym_mat3<int> a1 = a.tensor_transform(c),
                  correct_a1(198,18,1020,116,447,269);
    verify(__LINE__, a1, correct_a1);
    verify(__LINE__, a.tensor_transpose_transform(c),
                     sym_mat3<int>(371,9,967,148,597,238));
    sym_mat3<int> a1_bis;
    af::tiny<int, 6*6> m = c.tensor_transform_matrix();
    matrix_vector(6, 6, m.begin(), a.begin(), a1_bis.begin());
    verify(__LINE__, a1_bis, a1);
    vec3<int> v(3,-2,1);
    mat3<int> vc( 0, 3,-2,
                 -3, 0, 1,
                  2,-1, 0);
    verify(__LINE__, a.antisymmetric_tensor_transform(v),
                     sym_mat3<int>(-12,30,8,-9,10,-17));
    verify(__LINE__, a.tensor_transform(vc),
                     sym_mat3<int>(-12,30,8,-9,10,-17));
  }
  {
    mat3<int> a(3,4,-9, 2,1,7, -2,4,5);
    verify(__LINE__, a.self_times_self_transpose(),
                     sym_mat3<int>(106,54,45,-53,-35,35));
    verify(__LINE__, a.self_transpose_times_self(),
                     sym_mat3<int>(17,33,155,6,-23,-9));
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
