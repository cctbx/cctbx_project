#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/io.h>
#include <scitbx/sparse/lu_factorization.h>
#include <iostream>


void exercise_vector() {
  std::cout << "**** Exercising sparse::vector...\n\n";
  using namespace scitbx;
  typedef sparse::vector<double> sparse_vec;
  sparse_vec u(3);
    u[0] = 6.;
    u[2] = 8.;
  std::cout << sparse::dense_display(u);
  std::cout << std::endl;
  sparse_vec v(8);
    v[1] = 1.;
    v[4] = 2.;
    v[6] =-1.;
  std::cout << sparse::dense_display(v);
  std::cout << std::endl;
  std::cout << sparse::compressed_display(v);
  std::cout << std::endl;
  for (sparse_vec::iterator i=v.begin(); i != v.end(); i++) {
    std::cout << i.index() << ", ";
  }
  std::cout << std::endl;
  for (sparse_vec::iterator i=v.begin(); i != v.end() && i.index() < 5; i++) {
    std::cout << i.index() << ", ";
  }
  std::cout << std::endl << std::endl;
}

void exercise_matrix() {
  std::cout << "**** Exercising sparse::matrix...\n\n";
  using namespace scitbx;
  sparse::matrix<double> a(3,3);
  a(0,0) = 1.;
  a(0,2) = 2.;
  a(1,0) = 3.;
  a(2,1) = 4.;
  std::cout << "a = " << dense_display(a);
  sparse::vector<double> v1(3);
  v1[0] = -1.;
  std::cout << "\nv1 = " << dense_display(v1);
  std::cout << "\na v1 = " << dense_display(a*v1);
  sparse::vector<double> v2(3);
  v2[1] = 2.;
  v2[2] = 3.;
  std::cout << "\nv2 = " << dense_display(v2);
  std::cout << "\na v2 = " << dense_display(a*v2);
  std::cout << std::endl << std::endl;
}

static int lu_test_count = 0;

void assert_lu_factorization(const scitbx::sparse::matrix<double>& a) {
  std::cout << "**** Exercising sparse LU(" << ++lu_test_count << ")...\n\n";
  using namespace scitbx;
  sparse::gilbert_peierls_lu_factorization< sparse::matrix<double> > lu(a);
  std::cout << std::endl;
  std::cout << std::setw(6) << lu << "\n" << std::endl;
  sparse::matrix<double> a1 = lu.l() * lu.u(),
        a2 = lu.factored().deep_copy().permute_rows(lu.rows_permutation());
  SCITBX_ASSERT(lu.l().is_unit_lower_triangular());
  SCITBX_ASSERT(lu.u().is_upper_triangular());
  SCITBX_ASSERT(sparse::approx_equal(a1,a2))
               (sparse::dense_display(a1))(sparse::dense_display(a2));
  std::cout << "OK!\n" << std::endl;
}

void exercise_lu() {
  using namespace scitbx;

  /* mathematica
        a := SparseArray[ {
              {3, j_}  -> 1.5 - j/5,
              {7, j_}  -> -0.8 + j/5,
              {_, 5}   -> 2.1,
              {i_, i_} -> i
              }, {8, 8} ]
    */
  sparse::matrix<double> a(8,8);
  for (unsigned j=1; j <= a.n_cols(); j++) {
    sparse::vector<double> &c = a.col(j-1);
    for(unsigned i=1; i <= a.n_rows(); i++) {
      switch(i) {
        case 3:
          c[i-1] = 1.5 - j/5.;
          break;
        case 7:
          c[i-1] = -0.8 + j/5.;
          break;
        default:
          if      (j==5) c[i-1] = 2.1;
          else if (i==j) c[i-1] = i;
      }
    }
  }
  assert_lu_factorization(a);

  /* rectangular matrix m x n with m < n */
  sparse::matrix<double> b(5,8);
  b(4,0) = 1.;
  b(1,1) = -1.;
  b(1,2) = 0.5;
  b(2,1) = 1.8;
  b(2,2) = -2.;
  b(0,3) = 1.;
  b(2,4) = -1.;
  b(2,5) = 1.;
  b(2,6) = 0.5;
  b(3,5) = 0.5;
  b(3,6) = 1.;
  b(0,7) = 0.1;
  b(1,7) = 0.2;
  assert_lu_factorization(b);

  sparse::matrix<double> c = b.transpose();
  assert_lu_factorization(c);
}

int main() {
  exercise_vector();
  exercise_matrix();
  exercise_lu();
}
