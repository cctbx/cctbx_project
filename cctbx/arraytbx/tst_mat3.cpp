#include <cctbx/mat3.h>
#include <cctbx/constants.h>
#include <cctbx/array_family/simple_io.h>

using namespace cctbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int argc, char* argv[])
{
  {
    mat3<int> va;
    mat3<int> vb(0,0,0,0,0,0,0,0,0);
    mat3<int> vc(af::tiny_plain<int,9>(0,0,0,0,0,0,0,0,0));
    int id[] = {1,1,1,1,1,1,1,1,1};
    mat3<int> vd(id);
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
  }
  {
    mat3<double> r3_111(120. * constants::pi_180, vec3<double>(1.,1.,1.));
    mat3<double> r3_111i(-120. * constants::pi_180, vec3<double>(1.,1.,1.));
    {
      mat3<double> t = r3_111 * r3_111i - mat3<double>(1.);
      check_true(__LINE__, std::fabs(af::max(t.ref())) < 1.e-6);
    }
    {
      mat3<double> t = r3_111.inverse() - r3_111i;
      check_true(__LINE__, std::fabs(af::max(t.ref())) < 1.e-6);
    }
    {
      mat3<double> t = r3_111 - r3_111i.inverse();
      check_true(__LINE__, std::fabs(af::max(t.ref())) < 1.e-6);
    }
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

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
