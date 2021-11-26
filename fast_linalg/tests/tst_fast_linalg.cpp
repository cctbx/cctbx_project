#include <fast_linalg/cblas.h>
#include <fast_linalg/lapacke.h>
#include <fast_linalg/environment.h>
#include <iostream>

bool check_full_full(int N, const double *x, const double *y) {
  bool correct = true;
  for(int i=0; i<N; i++) for(int j=i; j<N; j++) {
    int k = j + N*i;
    if(x[k] != y[k]) {
      correct = false;
      std::cout << "\nExpected C[" << i << "," << j << "]=" << y[k] << "\n"
                << "but got C[" << i << "," << j << "]=" << x[k] << "\n";
    }
  }
  if(correct) std::cout << "passed!\n";
  return correct;
}

bool check_packed_packed(int N, const double *x, const double *y) {
  bool correct = true;
  for(int i=0; i<N*(N+1)/2; i++) {
    if(x[i] != y[i]) {
      correct = false;
      std::cout << "\nExpected " << y[i] << " at " << i << "-th position\n"
                << "but got " << x[i] << "\n";
    }
  }
  if(correct) std::cout << "passed!\n";
  return correct;
}

bool check_packed_full(int N, const double *x, const double *y) {
  bool correct = true;
  const double *px = x;
  for(int i=0; i<N; i++) for(int j=i; j<N; j++) {
    int k = j + N*i;
    if(y[k] != *px) {
      correct = false;
      std::cout << "\nExpected C[" << i << "," << j << "]=" << y[k] << "\n"
                << "but got C[" << i << "," << j << "]=" << *px << "\n";
    }
    px++;
  }
  if(correct) std::cout << "passed!\n";
  return correct;
}
using namespace fast_linalg;
int main() {
  bool all_correct = true;

  // First let's print some info about OpenBLAS
  fast_linalg::environment env;
  std::cout << "This library was build with the following options:\n";
  std::cout << env.build_config() << "\n";
  std::cout << "Running on " << env.physical_cores()
            << " " << env.cpu_family() << " cores\n";
  std::cout << "We will use " << env.threads() << " threads\n";

  const int N = 6;
  const int K = 4;
  double a[N*K] = { 1,  2,  3,  4,
                    5,  6,  7,  8,
                    9, 10, 11, 12,
                   13, 14, 15, 16,
                   17, 18, 19, 20,
                   21, 22, 23, 24};

  // C := A A^T with A being row-major and C being rectangular using CBLAS
  std::cout << "Rectangular Full Data # "
            << "[" << N << "x" << K << "]" << "[" << N << "x" << K << "]^T : ";
  double c[N*N];
  fast_linalg::syrk(
    CblasRowMajor, CblasUpper, CblasNoTrans, N, K,
              1.0, a, K, 0.0, c, N);
  double c_ref[N*N] = { 30,  70, 110, 150,  190,  230,
                         0, 174, 278, 382,  486,  590,
                         0,   0, 446, 614,  782,  950,
                         0,   0,   0, 846, 1078, 1310,
                         0,   0,   0,   0, 1374, 1670,
                         0,   0,   0,   0,    0, 2030};
  all_correct = all_correct && check_full_full(N, c, c_ref);

  // D := A A^T with A being column-major and D being RFP format using LAPACK
  std::cout << "Rectangular Full Packed # "
            << "[" << N << "x" << K << "]" << "[" << N << "x" << K << "]^T : ";
  // Usual trick: pretend we use column-major and demand transposed operation
  double cp[N*(N+1)/2];
  fast_linalg::sfrk(fast_linalg::LAPACK_COL_MAJOR, 'N', 'L', 'T', N, K, 1.0, a, K, 0.0, cp);
  /* Mathematica code to compute the expected result:
  c11 = ArrayFlatten[{{0}, { LowerTriangularize[c[[;; 3, ;; 3]]]}}] +
   ArrayFlatten[{{UpperTriangularize[c[[4 ;;, 4 ;;]]]}, {0}}];
  cpl = ArrayFlatten[{{c11}, {c[[4 ;;, ;; 3]]}}];
  cpl // MatrixForm
  Flatten[Transpose[cpl]]
  */
  double cp_ref[N*(N+1)/2] = { 846,   30,   70,
                               110,  150,  190,
                               230, 1078, 1374,
                               174,  278,  382,
                               486,  590, 1310,
                              1670, 2030,  446,
                               614,  782,  950};
  all_correct = all_correct && check_packed_packed(N, cp, cp_ref);

  // Copy cp into cp1 which is the usual packed upper triangle
  // The trick here is that the lower triangle in column-major
  // is the upper triangle in row-major that we want.
  double cp1[N*(N+1)/2];
  fast_linalg::tfttp(fast_linalg::LAPACK_COL_MAJOR, 'N', 'L', N, cp_ref, cp1);
  std::cout << "Rectangular full packed -> packed # ";
  all_correct = all_correct && check_packed_full(N, cp1, c_ref);

  if(all_correct) std::cout << "OK\n";
  return all_correct ? 0 : 1;
}
