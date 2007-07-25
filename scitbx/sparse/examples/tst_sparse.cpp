#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/io.h>
#include <scitbx/sparse/lu_factorization.h>
#include <iostream>


void exercise_vector() {
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
}

void exercise_matrix() {
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
}

void exercise_lu() {
  /* mathematica
    a := SparseArray[ { {3, j_}  -> 1.5 - j/5,
                        {7, j_}  -> -0.8 + j/5,
                        {_, 5}   -> 2.1,
                        {i_, i_} -> i },
                      {8, 8} ]
  */
  using namespace scitbx;
  std::cout << std::endl << std::endl;
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
  sparse::gilbert_peierls_lu_factorization< sparse::matrix<double> > lu(a);
  std::cout << std::setw(6) << lu;
  SCITBX_ASSERT(sparse::approx_equal(
                lu.l() * lu.u(),
                lu.factored().deep_copy().permute_rows(lu.rows_permutation()),
                1e-6));
}

int main() {
  exercise_vector();
  exercise_matrix();
  exercise_lu();
}
