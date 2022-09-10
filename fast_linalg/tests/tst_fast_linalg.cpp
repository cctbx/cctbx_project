#include <fast_linalg/cblas.h>
#include <fast_linalg/lapacke.h>
#include <fast_linalg/environment.h>
#include <iostream>
#include <stdio.h>

bool check_full_full(int N, const double* x, const double* y) {
  bool correct = true;
  for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
    int k = j + N * i;
    if (x[k] != y[k]) {
      correct = false;
      std::cout << "\nExpected C[" << i << "," << j << "]=" << y[k] << "\n"
        << "but got C[" << i << "," << j << "]=" << x[k] << "\n";
    }
  }
  if (correct) std::cout << "passed!\n";
  return correct;
}

bool check_packed_packed(int N, const double* x, const double* y) {
  bool correct = true;
  for (int i = 0; i < N * (N + 1) / 2; i++) {
    if (x[i] != y[i]) {
      correct = false;
      std::cout << "\nExpected " << y[i] << " at " << i << "-th position\n"
        << "but got " << x[i] << "\n";
    }
  }
  if (correct) std::cout << "passed!\n";
  return correct;
}

bool check_packed_full(int N, const double* x, const double* y) {
  bool correct = true;
  const double* px = x;
  for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
    int k = j + N * i;
    if (y[k] != *px) {
      correct = false;
      std::cout << "\nExpected C[" << i << "," << j << "]=" << y[k] << "\n"
        << "but got C[" << i << "," << j << "]=" << *px << "\n";
    }
    px++;
  }
  if (correct) std::cout << "passed!\n";
  return correct;
}
using namespace fast_linalg;
extern void print_matrix(char* desc, lapack_int m, lapack_int n, std::complex<double>* a, lapack_int lda);
void print_matrix(char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda);

template <typename FT>
bool cmp_arrays(const FT* a, const FT* b, size_t l) {
  int m = 100;
  for (size_t i = 0; i < l; i++) {
    if ((int)(a[i] * m) - (int)(b[i] * m) > 0) {
      printf("\n!!%6.2f != %6.2f", a[i], b[i]);
      return false;
    }
  }
  return true;
}

template <typename CT>
bool cmp_arraysc(const CT* a, const CT* b, size_t l) {
  int m = 100;
  for (size_t i = 0; i < l; i++) {
    if ((int)(a[i].real() * m) - (int)(b[i].real() * m) > 0) {
      printf("\n!!%6.2f != %6.2f", a[i].real(), b[i].real());
      return false;
    }
    if ((int)(a[i].imag() * m) - (int)(b[i].imag() * m) > 0) {
      printf("\n!!%6.2f != %6.2f", a[i].imag(), b[i].imag());
      return false;
    }
  }
  return true;
}

int main(int argc, char** argv) {
  bool all_correct = true;
  if (argc > 1) {
    printf("Initialising from %s\n", argv[1]);
    try {
      fast_linalg::initialise(argv[1]);
    }
    catch (...) {
    }
  }
  else {
    try {
      const char* ln = "libopenblas.dll";
      printf("Initialising from %s\n", ln);
      fast_linalg::initialise(ln);
    }
    catch (...) {
      try {
        const char* ln = "libopenblas.so";
        printf("Initialising from %s\n", ln);
        fast_linalg::initialise(ln);
      }
      catch (...) {
      }
    }
  }
  if (!fast_linalg::is_initialised()) {
    printf("Failed to initialise fast_linalg!\n");
    exit(1);
  }
  // First let's print some info about OpenBLAS
  fast_linalg::environment env;
  std::cout << "This library was build with the following options:\n";
  std::cout << env.build_config() << "\n";
  std::cout << "Running on " << env.physical_cores()
            << " " << env.cpu_family() << " cores\n";
  std::cout << "We will use " << env.threads() << " threads\n";
  {
    const int N = 6;
    const int K = 4;
    double a[N * K] = { 1,  2,  3,  4,
                      5,  6,  7,  8,
                      9, 10, 11, 12,
                     13, 14, 15, 16,
                     17, 18, 19, 20,
                     21, 22, 23, 24 };

    // C := A A^T with A being row-major and C being rectangular using CBLAS
    std::cout << "Rectangular Full Data # "
      << "[" << N << "x" << K << "]" << "[" << N << "x" << K << "]^T : ";
    double c[N * N];
    fast_linalg::syrk(
      CblasRowMajor, CblasUpper, CblasNoTrans, N, K,
      1.0, a, K, 0.0, c, N);
    double c_ref[N * N] = { 30,  70, 110, 150,  190,  230,
                           0, 174, 278, 382,  486,  590,
                           0,   0, 446, 614,  782,  950,
                           0,   0,   0, 846, 1078, 1310,
                           0,   0,   0,   0, 1374, 1670,
                           0,   0,   0,   0,    0, 2030 };
    all_correct = all_correct && check_full_full(N, c, c_ref);

    // D := A A^T with A being column-major and D being RFP format using LAPACK
    std::cout << "Rectangular Full Packed # "
      << "[" << N << "x" << K << "]" << "[" << N << "x" << K << "]^T : ";
    // Usual trick: pretend we use column-major and demand transposed operation
    double cp[N * (N + 1) / 2];
    fast_linalg::sfrk(fast_linalg::LAPACK_COL_MAJOR, 'N', 'L', 'T', N, K, 1.0, a, K, 0.0, cp);
    /* Mathematica code to compute the expected result:
    c11 = ArrayFlatten[{{0}, { LowerTriangularize[c[[;; 3, ;; 3]]]}}] +
     ArrayFlatten[{{UpperTriangularize[c[[4 ;;, 4 ;;]]]}, {0}}];
    cpl = ArrayFlatten[{{c11}, {c[[4 ;;, ;; 3]]}}];
    cpl // MatrixForm
    Flatten[Transpose[cpl]]
    */
    double cp_ref[N * (N + 1) / 2] = { 846,   30,   70,
                                 110,  150,  190,
                                 230, 1078, 1374,
                                 174,  278,  382,
                                 486,  590, 1310,
                                1670, 2030,  446,
                                 614,  782,  950 };
    if (!check_packed_packed(N, cp, cp_ref)) {
      all_correct = false;
    }

    // Copy cp into cp1 which is the usual packed upper triangle
    // The trick here is that the lower triangle in column-major
    // is the upper triangle in row-major that we want.
    double cp1[N * (N + 1) / 2];
    fast_linalg::tfttp(fast_linalg::LAPACK_COL_MAJOR, 'N', 'L', N, cp_ref, cp1);
    std::cout << "Rectangular full packed -> packed # ";
    if (!check_packed_full(N, cp1, c_ref)) {
      all_correct = false;
    }
  }
  {
    const lapack_int N = 4,
      LDA = N,
      LDVL = N,
      LDVR = N;
    lapack_int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR;
    /* Local arrays */
    std::complex<double> w[N], vl[LDVL * N], vr[LDVR * N];
    std::complex<double> a[LDA * N] = {
      {6.0161 , 1.5243}, {-3.5696 ,-6.9555}, {8.4695 ,-5.7720}, {7.8672 ,-6.2221},
      {3.5749 , 4.0301}, {-8.9314 , 7.1091}, {-1.6872 ,-6.7515}, {3.6823 , 7.5440},
      {1.6773 ,-3.7865}, {8.6006 , 5.4106}, {-1.6935 ,-9.8013}, {-6.2401 ,-6.4186},
      {6.7681 ,-4.5137}, {-5.6985 , 8.0920}, {2.1531 ,-4.8516}, {-8.7264 ,-0.1071}
    };
    printf("\nTest geev - general matrix eigen problem.\n");
    lapack_int info = fast_linalg::geev(LAPACK_ROW_MAJOR, LAPACK_EIGENVALUES_AND_EIGENVECTORS,
      LAPACK_EIGENVALUES_AND_EIGENVECTORS, n, a, lda, w, vl,
      ldvl, vr, ldvr);
    if (info > 0) {
      printf("Failed to compute eigenvalues.\n");
      all_correct = false;
    }
    else {
      print_matrix("Eigenvalues", 1, n, w, 1);
      print_matrix("Left eigenvectors", n, n, vl, ldvl);
      print_matrix("Right eigenvectors", n, n, vr, ldvr);
      std::complex<double> w_tst[] = {
        {2.3797,-14.9575}, {8.3883,-2.7927}, {-11.9239,1.2782}, {-12.1793,15.1970} };
      std::complex<double> vl_tst[] = {
         {0.3059,0.0463}, {0.6846,0.0000}, {-0.1280,-0.3516}, {-0.1127,0.1348},
         {-0.0363,-0.3859}, {0.1345,0.3407}, {-0.4961,0.0883}, {0.8112,0.0000},
        {0.7488,0.0000}, {-0.0083,0.5602}, {0.0598,0.3013}, {-0.2160,-0.2566},
        {0.0752,0.4332}, {0.2884,0.0103}, {0.7154,0.0000}, {0.2851,0.3425} };
      std::complex<double> vr_tst[] = {
        {0.0300,0.5231}, {0.8717,0.0000}, {0.6162,0.0000}, {-0.2737,0.1951},
        {-0.0563,-0.1815}, {0.1617,0.3467}, {0.1741,-0.1522}, {0.6838,0.0000},
        {0.7645,0.0000}, {-0.0686,0.1140}, {-0.1895,-0.3938}, {-0.0635,-0.1354},
        {0.1819,0.2680}, {0.1830,-0.2064}, {-0.6003,-0.1244}, {0.5834,0.2383} };
      if (all_correct) {
        all_correct = cmp_arraysc(&w_tst[0], &w[0], N);
      }
      if (all_correct) {
        all_correct = cmp_arraysc(&vl_tst[0], &vl[0], LDVL * N);
      }
      if (all_correct) {
        all_correct = cmp_arraysc(&vr_tst[0], &vr[0], LDVR * N);
      }
    }
    printf("\nTest getrf/getri - general matrix eigen problem.\n");
    lapack_int ipiv[N];
    info = fast_linalg::getrf(LAPACK_ROW_MAJOR, N, N, a, LDA, ipiv);
    if (info > 0) {
      printf("Failed getrf - LU decomposition.\n");
      all_correct = false;
    }
    else {
      print_matrix("After getrf", N, N, a, 1);
    }
    info = fast_linalg::getri(LAPACK_ROW_MAJOR, N, a, LDA, ipiv);
    if (info > 0) {
      printf("Failed getri - matrix inverse.\n");
      all_correct = false;
    }
    else {
      print_matrix("After getrf", N, N, a, 1);
    }
  }
  {
    const lapack_int N = 4,
      LDA = N;
    lapack_int n = N, lda = LDA;
    double w[N];
    std::complex<double> a[LDA * N] = {
      {-8.1415, 0} ,      {0,0} ,            {0,0},             {0,0},
      {-1.5971,1.3513} ,  {0.6667,0} ,       {0,0} ,            {0,0},
      {-0.8081, 0.5664},  {0.5465, -1.8765}, {8.6014,0},        {0,0},
      {5.6678,-2.4521} ,  {-6.7320, 9.8121}, {1.7920, -5.4521}, {7.6382,0}
    };
    printf("\nTest heev - Hermitian matrix eigen problem\n");
    lapack_int info = fast_linalg::heev(LAPACK_ROW_MAJOR,
      LAPACK_EIGENVALUES_AND_EIGENVECTORS, LAPACK_LOWER, n, a, lda, w);
    if (info > 0) {
      printf("Failed to compute eigenvalues.\n");
      all_correct = false;
    }
    else {
      print_matrix("Eigenvalues", 1, n, w, 1);
      print_matrix("Eigenvectors (stored columnwise)", n, n, a, lda);

      double w_tst[] = { -12.8458, -5.5622, 7.3574, 19.8155 };
      std::complex<double> v_tst[] = {
        {-0.7032,0.0000}, {-0.6776,0.0000}, {0.1122,0.0000}, {-0.1840,0.0000},
        {0.2135,0.4198}, {-0.4081,-0.4845}, {-0.3705,-0.0302}, {0.4611,0.1616},
        {-0.1296,-0.0949}, {0.0602,0.1302}, {-0.8117,-0.3629}, {-0.2210,-0.3380},
        {0.5024,-0.0747}, {-0.3393,-0.0625}, {0.0258,-0.2409}, {-0.6547,0.3684} };
      if (all_correct) {
        all_correct = cmp_arrays(&w_tst[0], &w[0], N);
      }
      if (all_correct) {
        all_correct = cmp_arraysc(&v_tst[0], &a[0], N*N);
      }
    }
  }
  {
    const lapack_int N = 4,
      LDA = N;
    lapack_int n = N, lda = LDA;
    double w[N];
    double a[LDA * N] = {
      -8.7375, -6.2604, -0.8461,  9.1861,
      0,       7.6420,  -8.9406, -5.9790,
      0,       0,       2.5771,  3.8387,
      0,       0,       0,     - 2.7090};
    printf("\nTest syev - real symmetric eigen problem\n");
    lapack_int info = fast_linalg::syev(LAPACK_ROW_MAJOR,
      LAPACK_EIGENVALUES_AND_EIGENVECTORS, LAPACK_UPPER, n, a, lda, w);
    if (info > 0) {
      printf("Failed to compute eigenvalues.\n");
      all_correct = false;
    }
    else {
      print_matrix("Eigenvalues", 1, n, w, 1);
      print_matrix("Eigenvectors (columnwise)", n, n, a, lda);
      double w_tst[] = { -16.3168, -4.9809, 0.9332, 19.1371 };
      double v_tst[] = {
        0.8015, -0.0811, -0.5208, 0.2823,
        0.1581, -0.6618, -0.0499, -0.7311,
        0.2191, -0.5216, 0.6751, 0.4735,
        -0.5334, -0.5323, -0.5201, 0.4020 };
      if (all_correct) {
        all_correct = cmp_arrays(&w_tst[0], &w[0], N);
      }
      if (all_correct) {
        all_correct = cmp_arrays(&v_tst[0], &a[0], N*N);
      }
    }
  }
  if (all_correct) {
    std::cout << "OK\n";
  }
  return all_correct ? 0 : 1;
}

void print_matrix(char* desc, lapack_int m, lapack_int n, std::complex<double>* a, lapack_int lda) {
  printf("\n %s\n", desc);
  for (lapack_int i = 0; i < m; i++) {
    for (lapack_int j = 0; j < n; j++) {
      printf(" {%6.4f,%6.4f},", a[i * lda + j].real(), a[i * lda + j].imag());
    }
    printf("\n");
  }
}

void print_matrix(char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda) {
  printf("\n %s\n", desc);
  for (lapack_int i = 0; i < m; i++) {
    for (lapack_int j = 0; j < n; j++) {
      printf(" %6.4f,", a[i * lda + j]);
    }
    printf("\n");
  }
}