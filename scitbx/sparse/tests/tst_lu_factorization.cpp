#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/io.h>
#include <scitbx/sparse/lu_factorization.h>
#include <iostream>


static int lu_test_count = 0;

void exercise_lu_factorization(const scitbx::sparse::matrix<double>& a) {
  using namespace scitbx;
  sparse::gilbert_peierls_lu_factorization< sparse::matrix<double> > lu(a);
  bool failed = false;
  if (failed = failed || !lu.l().is_unit_lower_triangular())
    std::cout << "Error: L is not unit lower triangular\n";
  if (failed = failed || !lu.u().is_upper_triangular())
    std::cout << "Error: U is not upper triangular\n";
  sparse::matrix<double>
    a1 = lu.l() * lu.u(),
    a2 = lu.factored().deep_copy().permute_rows(lu.rows_permutation());
  if (failed = failed || !sparse::approx_equal(a1,a2))
    std::cout << "Error: PA != LU\n" << std::endl;
  if (failed) {
    std::cout << std::setw(6) << lu << "\n" << std::endl;
    std::cout << "PA = " << sparse::dense_display(a1) << "\n" << std::endl;
    std::cout << "LU = " << sparse::dense_display(a2) << "\n" << std::endl;
  }
  else std::cout << "sparse LU(" << ++lu_test_count << "): OK\n";
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
  exercise_lu_factorization(a);

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
  exercise_lu_factorization(b);

  sparse::matrix<double> c = b.transpose();
  exercise_lu_factorization(c);
}

int main() {
  exercise_lu();
  return 0;
}
