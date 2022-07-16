#include <fast_linalg/cblas.h>
#include <fast_linalg/lapacke.h>
#include <fast_linalg/environment.h>
#include <iostream>

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
extern void print_matrix(char* desc, lapack_int m, lapack_int n, _lapack_complex_double* a, lapack_int lda);
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
    if ((int)(a[i].real * m) - (int)(b[i].real * m) > 0) {
      printf("\n!!%6.2f != %6.2f", a[i].real, b[i].real);
      return false;
    }
    if ((int)(a[i].imag * m) - (int)(b[i].imag * m) > 0) {
      printf("\n!!%6.2f != %6.2f", a[i].imag, b[i].imag);
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
    all_correct = all_correct && check_packed_packed(N, cp, cp_ref);

    // Copy cp into cp1 which is the usual packed upper triangle
    // The trick here is that the lower triangle in column-major
    // is the upper triangle in row-major that we want.
    double cp1[N * (N + 1) / 2];
    fast_linalg::tfttp(fast_linalg::LAPACK_COL_MAJOR, 'N', 'L', N, cp_ref, cp1);
    std::cout << "Rectangular full packed -> packed # ";
    all_correct = all_correct && check_packed_full(N, cp1, c_ref);
  }
  //https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/nonsymmetric-eigenproblems/geev-function/zgeev-example/lapacke-zgeev-example-c-column.html
  {
    const lapack_int N = 4,
      LDA = N,
      LDVL = N,
      LDVR = N;
    lapack_int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR;
    /* Local arrays */
    _lapack_complex_double w[N], vl[LDVL * N], vr[LDVR * N];
    _lapack_complex_double a[LDA * N] = {
           {-3.84,  2.25}, {-8.94, -4.75}, { 8.95, -6.53}, {-9.87,  4.82},
           {-0.66,  0.83}, {-4.40, -3.82}, {-3.50, -4.26}, {-3.15,  7.36},
           {-3.99, -4.73}, {-5.88, -6.60}, {-3.36, -0.40}, {-0.75,  5.23},
           { 7.74,  4.18}, { 3.66, -7.53}, { 2.58,  3.60}, { 4.59,  5.41}
    };
    printf("\nLAPACKE_zgeev. Example Program Results\n");
    lapack_int info = fast_linalg::geev(LAPACK_ROW_MAJOR, LAPACK_EIGENVALUES_AND_EIGENVECTORS,
      LAPACK_EIGENVALUES_AND_EIGENVECTORS, n, a, lda, w, vl,
      ldvl, vr, ldvr);
    if (info > 0) {
      printf("The algorithm failed to compute eigenvalues.\n");
      all_correct = false;
    }
    else {
      print_matrix("Eigenvalues", 1, n, w, 1);
      print_matrix("Left eigenvectors", n, n, vl, ldvl);
      print_matrix("Right eigenvectors", n, n, vr, ldvr);
      _lapack_complex_double w_tst[] = {
        {-9.4299,-12.9833}, {-3.4418,12.6897}, {0.1055,-3.3950}, {5.7562,7.1286} };
      _lapack_complex_double vl_tst[] = {
        {0.2414,-0.1847}, {0.6135,0.0000}, {-0.1828,-0.3347}, {0.2765,0.0884},
        {0.7861,0.0000}, {-0.0499,-0.2721}, {0.8218,0.0000}, {-0.5477,0.1572},
        {0.2195,-0.2689}, {-0.2088,0.5347}, {-0.3714,0.1525}, {0.4451,0.0912},
        {-0.0170,0.4109}, {0.4027,-0.2353}, {0.0575,0.1208}, {0.6202,0.0000} };
      _lapack_complex_double vr_tst[] = {
        {0.4309,0.3268}, {0.8257,0.0000}, {0.5984,0.0000}, {-0.3054,0.0333},
        {0.5087,-0.0288}, {0.0750,-0.2487}, {-0.4005,-0.2014}, {0.0398,0.3445},
        {0.6198,0.0000}, {-0.2458,0.2789}, {-0.0901,-0.4753}, {0.3583,0.0606},
        {-0.2269,0.1104}, {-0.1034,-0.3192}, {-0.4348,0.1337}, {0.8082,0.0000} };
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
  }
  //https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/symmetric-eigenproblems/heev-function/zheev-example/lapacke-zheev-example-c-row.html
  {
    const lapack_int N = 4,
      LDA = N;
    lapack_int n = N, lda = LDA;
    double w[N];
    _lapack_complex_double a[LDA * N] = {
       { 9.14,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
       {-4.37,  9.22}, {-3.35,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
       {-1.98,  1.72}, { 2.25,  9.51}, {-4.82,  0.00}, { 0.00,  0.00},
       {-8.96,  9.50}, { 2.57, -2.40}, {-3.24, -2.04}, { 8.44,  0.00}
    };
    printf("\nLAPACKE_zheev (row-major, high-level) Example Program Results\n");
    lapack_int info = fast_linalg::heev(LAPACK_ROW_MAJOR,
      LAPACK_EIGENVALUES_AND_EIGENVECTORS, 'L', n, a, lda, w);
    if (info > 0) {
      printf("The algorithm failed to compute eigenvalues.\n");
      all_correct = false;
    }
    else {
      print_matrix("Eigenvalues", 1, n, w, 1);
      print_matrix("Eigenvectors (stored columnwise)", n, n, a, lda);

      double w_tst[] = { -16.0047, -6.7650, 6.6657, 25.5140 };
      _lapack_complex_double v_tst[] = {
        {0.3448,-0.0000}, {-0.5457,0.0000}, {0.3112,0.0000}, {-0.6975,0.0000},
        {0.4418,-0.5389}, {0.2620,0.1810}, {0.4536,0.2868}, {0.2158,-0.2800},
        {-0.4795,-0.3744}, {-0.5195,-0.0157}, {-0.0524,0.5734}, {0.1461,0.0830},
        {0.1005,-0.1236}, {-0.5030,0.2787}, {-0.2285,-0.4810}, {0.3413,-0.4938} };
      if (all_correct) {
        all_correct = cmp_arrays(&w_tst[0], &w[0], N);
      }
      if (all_correct) {
        all_correct = cmp_arraysc(&v_tst[0], &a[0], N*N);
      }
    }
  }
  {
    const lapack_int N = 5,
      LDA = N;
    lapack_int n = N, lda = LDA;
    double w[N];
    double a[LDA * N] = {
        1.96, -6.49, -0.47, -7.20, -0.65,
        0.00,  3.80, -6.39,  1.50, -6.34,
        0.00,  0.00, 4.17, -1.51, 2.67,
        0.00,  0.00, 0.00,  5.70, 1.80,
        0.00,  0.00, 0.00,  0.00, -7.10
    };
    printf("\nLAPACKE_dsyev (row-major, high-level) Example Program Results\n");
    lapack_int info = fast_linalg::syev(LAPACK_ROW_MAJOR,
      LAPACK_EIGENVALUES_AND_EIGENVECTORS, 'U', n, a, lda, w);
    if (info > 0) {
      printf("The algorithm failed to compute eigenvalues.\n");
      all_correct = false;
    }
    else {
      print_matrix("Eigenvalues", 1, n, w, 1);
      print_matrix("Eigenvectors (stored columnwise)", n, n, a, lda);
      double w_tst[] = { -11.0656, -6.2287, 0.8640, 8.8655, 16.0948 };
      double v_tst[] = {
        -0.2981, -0.6075, 0.4026, -0.3745, 0.4896,
        -0.5078, -0.2880, -0.4066, -0.3572, -0.6053,
        -0.0816, -0.3843, -0.6600, 0.5008, 0.3991,
        -0.0036, -0.4467, 0.4553, 0.6204, -0.4564,
        -0.8041, 0.4480, 0.1725, 0.3108, 0.1622 };
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

void print_matrix(char* desc, lapack_int m, lapack_int n, _lapack_complex_double* a, lapack_int lda) {
  printf("\n %s\n", desc);
  for (lapack_int i = 0; i < m; i++) {
    for (lapack_int j = 0; j < n; j++) {
      printf(" {%6.4f,%6.4f},", a[i * lda + j].real, a[i * lda + j].imag);
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