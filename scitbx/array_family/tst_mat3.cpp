#include <cmath>
#include <scitbx/mat3.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/simple_io.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int /*argc*/, char* /*argv*/[])
{
  {
    mat3<int> va;
    check_true(__LINE__, va.size() == 9);
    mat3<int> vb(0,0,0,0,0,0,0,0,0);
    for(int j=0;j<9;j++) check_true(__LINE__, vb[j] == 0);
    mat3<int> vc(af::tiny_plain<int,9>(0,0,0,0,0,0,0,0,0));
    for(int j=0;j<9;j++) check_true(__LINE__, vc[j] == 0);
    int id[] = {1,1,1,1,1,1,1,1,1};
    mat3<int> vd(id);
    for(int j=0;j<9;j++) check_true(__LINE__, vd[j] == 1);
  }
  {
    mat3<double> diag9(1,0,0,0,1,0,0,0,1);
    mat3<double> diag1(1.0);
    mat3<double> diag3(1,1,1);
    check_true(__LINE__, diag1 == diag9);
    check_true(__LINE__, diag3 == diag9);
    mat3<double> diag234(2,3,-4);
    mat3<double> diag2349(2,0,0,0,3,0,0,0,-4);
    check_true(__LINE__, diag234 == diag2349);
  }
  {
    mat3<int> a; a.fill(0);
    mat3<int> b; b.fill(0);
    mat3<int> c; c.fill(1);
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
    check_true(__LINE__, mat3<int>(1) == mat3<int>(1,0,0,0,1,0,0,0,1));
    check_true(__LINE__, mat3<int>(vec3<int>(1,2,3)) == mat3<int>(
      1,0,0,0,2,0,0,0,3));
  }
  {
    mat3<int> a(1,2,3,4,5,6,7,8,9);
    mat3<int> b(4,5,6,7,8,9,10,11,12);
    mat3<int> c(5,7,9,11,13,15,17,19,21);
    mat3<int> d(2,4,6,8,10,12,14,16,18);
    mat3<int> r3(0,-1,0, 1,-1,0, 0,0,1);
    mat3<int> r3i(-1,1,0, -1,0,0, 0,0,1);
    check_true(__LINE__, a(0,0) == 1);
    check_true(__LINE__, a(0,1) == 2);
    check_true(__LINE__, a(0,2) == 3);
    check_true(__LINE__, a(1,0) == 4);
    check_true(__LINE__, a(1,1) == 5);
    check_true(__LINE__, a(1,2) == 6);
    check_true(__LINE__, a(2,0) == 7);
    check_true(__LINE__, a(2,1) == 8);
    check_true(__LINE__, a(2,2) == 9);
    check_true(__LINE__, a + b == c);
    check_true(__LINE__, a + 3 == b);
    check_true(__LINE__, 3 + a == b);
    check_true(__LINE__, a - b == -3);
    check_true(__LINE__, b - 3 == a);
    check_true(__LINE__, c - b == a);
    check_true(__LINE__, a * b == mat3<int>(
      48,54,60,111,126,141,174,198,222));
    check_true(__LINE__, transpose_mul(a, b) == mat3<int>(
      102,114,126,123,138,153,144,162,180));
    check_true(__LINE__, mul_transpose(a, b) == mat3<int>(
      32,50,68,77,122,167,122,194,266));
    check_true(__LINE__, transpose_mul(a.transpose(), b) == mat3<int>(
      48,54,60,111,126,141,174,198,222));
    check_true(__LINE__, mul_transpose(a, b.transpose()) == mat3<int>(
      48,54,60,111,126,141,174,198,222));
    check_true(__LINE__, transpose_mul(a, b.transpose()) == mat3<int>(
      66,102,138,81,126,171,96,150,204));
    check_true(__LINE__, mul_transpose(a.transpose(), b) == mat3<int>(
      66,102,138,81,126,171,96,150,204));
    check_true(__LINE__, a * vec3<int>(1,2,3) == vec3<int>(
      14,32,50));
    check_true(__LINE__, std::fabs(af::max(
      (a * vec3<double>(1,2,3) - vec3<double>(14,32,50)).ref())) < 1.e-6);
    check_true(__LINE__, vec3<int>(1,2,3) * a == vec3<int>(
      30,36,42));
    check_true(__LINE__, std::fabs(af::max(
      (vec3<double>(1,2,3) * a - vec3<double>(30,36,42)).ref())) < 1.e-6);
    check_true(__LINE__, a * 2 == d);
    check_true(__LINE__, 2 * a == d);
    check_true(__LINE__, d / 2 == a);
    check_true(__LINE__, 5040 / d == mat3<int>(
      2520,1260,840,630,504,420,360,315,280));
    check_true(__LINE__, a % 2 == mat3<int>(1,0,1,0,1,0,1,0,1));
    check_true(__LINE__, 4 % d == mat3<int>(0,0,4,4,4,4,4,4,4));
    {
      mat3<int> t(a);
      check_true(__LINE__, (t += b) == c);
    }
    {
      mat3<int> t(a);
      check_true(__LINE__, (t += 3) == b);
    }
    {
      mat3<int> t(a);
      check_true(__LINE__, (t -= b) == -3);
    }
    {
      mat3<int> t(b);
      check_true(__LINE__, (t -= 3) == a);
    }
    {
      mat3<int> t(a);
      check_true(__LINE__, (t *= 2) == d);
    }
    {
      mat3<int> t(d);
      check_true(__LINE__, (t /= 2) == a);
    }
    {
      mat3<int> t(a);
      check_true(__LINE__, (t %= 2) == mat3<int>(1,0,1,0,1,0,1,0,1));
    }
    check_true(__LINE__, -a == mat3<int>(-1,-2,-3,-4,-5,-6,-7,-8,-9));
    check_true(__LINE__, +a == a);
    check_true(__LINE__, a.get_row(0) == vec3<int>(1,2,3));
    check_true(__LINE__, a.get_row(1) == vec3<int>(4,5,6));
    check_true(__LINE__, a.get_row(2) == vec3<int>(7,8,9));
    {
      mat3<int> t; t.fill(0);
      t.set_row(0, vec3<int>(1,2,3));
      t.set_row(1, vec3<int>(4,5,6));
      t.set_row(2, vec3<int>(7,8,9));
      check_true(__LINE__, t == a);
      t.set_row(0, vec3<int>(4,5,6));
      t.set_row(1, vec3<int>(1,2,3));
      t.swap_rows(0, 1);
      check_true(__LINE__, t == a);
      t.set_row(1, vec3<int>(7,8,9));
      t.set_row(2, vec3<int>(4,5,6));
      t.swap_rows(1, 2);
      check_true(__LINE__, t == a);
    }
    check_true(__LINE__, a.get_column(0) == vec3<int>(1,4,7));
    check_true(__LINE__, a.get_column(1) == vec3<int>(2,5,8));
    check_true(__LINE__, a.get_column(2) == vec3<int>(3,6,9));
    {
      mat3<int> t; t.fill(0);
      t.set_column(0, vec3<int>(1,4,7));
      t.set_column(1, vec3<int>(2,5,8));
      t.set_column(2, vec3<int>(3,6,9));
      check_true(__LINE__, t == a);
      t.set_column(0, vec3<int>(2,5,8));
      t.set_column(1, vec3<int>(1,4,7));
      t.swap_columns(0, 1);
      check_true(__LINE__, t == a);
      t.set_column(1, vec3<int>(3,6,9));
      t.set_column(2, vec3<int>(2,5,8));
      t.swap_columns(1, 2);
      check_true(__LINE__, t == a);
    }
    check_true(__LINE__, a.diagonal() == vec3<int>(1,5,9));
    check_true(__LINE__, a.transpose() == mat3<int>(1,4,7,2,5,8,3,6,9));
    check_true(__LINE__, a.trace() == 15);
    check_true(__LINE__, a.determinant() == 0);
    check_true(__LINE__, mat3<int>(1).determinant() == 1);
    check_true(__LINE__, mat3<int>(2).determinant() == 8);
    check_true(__LINE__, r3.determinant() == 1);
    check_true(__LINE__, r3.co_factor_matrix_transposed() == r3i);
    check_true(__LINE__, r3.inverse() == r3i);
    check_true(__LINE__, a.dot(b) == 420);
    check_true(__LINE__, b.dot(a) == 420);
  }
  {
    check_true(__LINE__, mat3<int>(7,-2,9,-4,5,-6,1,-8,3).max_abs() == 9);
    check_true(__LINE__, mat3<int>(-7,2,-9,4,-5,6,-1,8,-3).max_abs() == 9);
  }
  {
    mat3<double> t;
    t.set_row(0, vec3<int>(1,2,3));
    t.set_row(1, vec3<int>(2,5,6));
    t.set_row(2, vec3<int>(3,6,9));
    check_true(__LINE__, t.is_symmetric());
    t[1] = 2.01;
    check_false(__LINE__, t.is_symmetric());
    check_true(__LINE__, t.is_symmetric(0.1));
    check_false(__LINE__, t.is_symmetric(0.001));
  }
  {
    mat3<int> a;
    a.set_row(0, vec3<int>(1,0,0));
    a.set_row(1, vec3<int>(0,5,0));
    a.set_row(2, vec3<int>(0,0,9));
    check_true(__LINE__, a.is_diagonal());
    a[1] = 2;
    check_false(__LINE__, a.is_diagonal());
  }
  {
    mat3<int> a(1,2,3,4,5,6,7,8,9);
    a.scale(vec3<int>(1,2,3));
    check_true(__LINE__, a == mat3<int>(1,4,9,4,10,18,7,16,27));
  }
  {
    mat3<double> a(3,2,3,5,7,3,4,2,7);
    check_true(__LINE__, std::fabs(af::max((a.ortho() - mat3<double>(
      3., -0.94, -1.16258352,
      5., 2.1, 0.12917595,
      4., -1.92, 0.71046771)).ref())) < 1.e-6);
    std::pair<mat3<double>, vec3<double> > d = a.decompose();
    check_true(__LINE__, std::fabs(af::max((d.first - mat3<double>(
      0.424264069, -0.313682063, -0.849472521,
      0.707106781,  0.700779076,  0.094385836,
      0.565685425, -0.640712298,  0.519122096)).ref())) < 1.e-6);
    check_true(__LINE__, std::fabs(af::max((d.second - vec3<double>(
      7.071067812,  2.996664813,  1.368594617)).ref())) < 1.e-6);
  }
  {
    vec3<int> a(2,5,11);
    vec3<int> b(-3,13,17);
    mat3<int> ax(mat3<int>::cross_product_matrix(a));
    mat3<int> bx(mat3<int>::cross_product_matrix(b));
    check_true(__LINE__, ax * b == a.cross(b));
    check_true(__LINE__, bx * a == b.cross(a));
    check_true(__LINE__, ax * a == a.cross(a));
    check_true(__LINE__, bx * b == b.cross(b));
    check_true(__LINE__, ax * a != a.cross(b));
    check_true(__LINE__, bx * b != b.cross(a));
    check_true(__LINE__, ax * b != a.cross(a));
    check_true(__LINE__, bx * a != b.cross(b));
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
