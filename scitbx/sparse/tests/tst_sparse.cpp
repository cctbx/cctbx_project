#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/array_family/initialiser.h>
#include <iostream>

namespace scitbx { namespace sparse {

  /// Test the efficient assignment and augmented assignment available in C++
  void exercise_vector_element_assignment() {
    {
      unsigned n = 4;
      vector<double> v(n);
      SCITBX_ASSERT(v.size() == n)(v.size());
      SCITBX_ASSERT(!v.is_compact());
      for (int i=0; i<n; ++i) SCITBX_ASSERT(v[i] == 0)(v[i]);
      v.compact();
      SCITBX_ASSERT(v.is_compact());
    }
    {
      unsigned n = 5;
      vector<double> v(n);
      v[1] += 1e-20;
      v[2] += 2.;
      v.compact();
      af::shared<double> w(n);
      w[2] = 2; w[1] = 1e-20;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
    }
    {
      unsigned n = 6;
      vector<double> v(n);
      v[4] += 1.;
      v[5] += 2.;
      v[4] += 3.;
      v[5] = 4.;
      v[3] = 5.;
      v[5] += 6.;
      SCITBX_ASSERT(v.size() == n);
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == n);
      af::shared<double> w(n);
      w[4] = 4.; w[5] = 10.; w[3] = 5.;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
      v[1] = 0;
      SCITBX_ASSERT(!v.is_compact());
    }
    {
      unsigned n = 10;
      vector<double> v(n);
      v[2] += 1.;
      v[3] += 2.;
      v[5] += 3.;
      v[3] += 4.;
      v[2] += 5.;
      v[5] += 6.;
      SCITBX_ASSERT(v.size() == n);
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == n);
      af::shared<double> w(n);
      w[2] = 6.; w[3] = 6.; w[5] = 9.;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
    }
    {
      unsigned n = 8;
      vector<double> v(n);
      SCITBX_ASSERT(v.size() == n);
      v[1] = 1.;
      v[4] = 2.;
      v[6] = 3.;
      SCITBX_ASSERT(v.size() == n);
      SCITBX_ASSERT(!v.is_compact());
      v.compact();
      SCITBX_ASSERT(v.is_compact());
      SCITBX_ASSERT(v.size() == n);
      af::shared<double> w(n);
      w[1] = 1.; w[4] = 2.; w[6] = 3.;
      for (vector<double>::iterator p=v.begin(); p != v.end(); ++p) {
        SCITBX_ASSERT(*p == w[p.index()])(p.index())(*p);
      }
    }
  }

  /// Permuted vector expression
  void exercise_permutation() {
    unsigned const n = 5;
    vector<double> u(n);
    u[1] = 1.;
    u[4] = 4.;
    std::vector<int> p(n);
    af::init(p) = 3, 2, 4, 1, 0;
    vector<double> v = permute(u, p);
    vector<double> v0(n);
    v0[2] = 1.;
    v0[0] = 4.;
    SCITBX_ASSERT(v == v0);
    u.permute(p);
    SCITBX_ASSERT(u == v);
  }


  /// Test vector linear combinations
  void exercise_vector_linear_algebra() {
    {
      unsigned n = 5;
      vector<double> u(n), v(n);
      u[0] = 1.;
      u[3] = 2.;
      v[1] = -1.;
      v[3] = 3.;
      v[4] = -2.;
      vector<double> w1 = 2.*u - v;
      vector<double> ww(n);
      ww[0] = 2.;
      ww[1] = 1.;
      ww[3] = 6.;
      ww[4] = 2.;
      SCITBX_ASSERT(w1 == ww);
      vector<double> w2 = -v + 2.*u;
      SCITBX_ASSERT(w2 == ww);
      vector<double> w3 = u+v;
      vector<double> ww3(n);
      ww3[0] = 1.;
      ww3[1] = -1.;
      ww3[3] = 5.;
      ww3[4] = -2.;
      SCITBX_ASSERT(w3 == ww3);
      vector<double> w4 = -2.*u + 3.*v;
      vector<double> ww4(n);
      ww4[0] = -2.;
      ww4[1] = -3.;
      ww4[3] = 5.;
      ww4[4] = -6.;
    }
  }

  void exercise_matrix_times_dense_vector() {
    int const m=5, n=3;
    double x_[n] = { 1., 2., 3. };
    double y_[m];
    af::ref<double> x(x_, n), y(y_, m);
    matrix<double> a(m, n);
    a(0, 0) =  1.;
    a(3, 0) =  -1.;
    a(1, 1) = -1.;
    a(2, 1) = -1.;
    a(0, 2) =  1.;
    a(4, 2) = -1.;
    a.set_compact(true);
    y = a*x;
    double y0_[m] = { 4, -2, -2, -1, -3 };
    af::ref<double> y0(y0_, m);
    SCITBX_ASSERT(y.all_eq(y0));
  }

  void exercise_matrix_transpose_times_dense_vector() {
    int const m=5, n=3;
    double x_[m] = { 1., 2., 3., 4., 5. };
    double y_[n];
    af::ref<double> x(x_, m), y(y_, n);
    matrix<double> a(m, n);
    a(0, 0) = -1.;
    a(3, 0) =  1.;
    a(1, 1) = -1.;
    a(2, 1) =  1.;
    a(0, 2) = -1.;
    a(4, 2) =  1.;
    a.set_compact(true);
    y = a.transpose_times(x);
    double y0_[n] = { 3, 1, 4 };
    af::ref<double> y0(y0_, n);
    SCITBX_ASSERT(y.all_eq(y0));
  }

  void exercise_operations_with_dense_matrices() {
    int const m=5, n=3;
    matrix<double> a(m,n);
    a(1, 0) = -1;
    a(2, 1) =  1;
    a(3, 2) =  2;
    af::versa<double, af::mat_grid> c0(af::mat_grid(m, n));
    c0(1, 0) = -1;
    c0(2, 1) =  1;
    c0(3, 2) =  2;
    af::versa<double, af::mat_grid> c1 = c0.deep_copy();
    c1(4, 1) = 5;
    af::versa<double, af::mat_grid> b0 = c1.deep_copy();
    b0 = a;
    SCITBX_ASSERT(b0.all_eq(c0));
    af::versa<double, af::mat_grid> d0 = a;
    SCITBX_ASSERT(d0.all_eq(c0));
    af::versa<double, af::mat_grid> b1(af::mat_grid(m, n));
    b1(4, 1) = 5;
    b1 += a;
    SCITBX_ASSERT(b1.all_eq(c1));
  }

}}

int main() {
  using namespace scitbx::sparse;
  exercise_operations_with_dense_matrices();
  exercise_vector_element_assignment();
  exercise_permutation();
  exercise_vector_linear_algebra();
  exercise_matrix_times_dense_vector();
  exercise_matrix_transpose_times_dense_vector();
  std::cout << "OK\n" << std::endl;
  return 0;
}
