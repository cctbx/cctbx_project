#include <fast_linalg/cblas.h>
#include <fast_linalg/lapacke.h>
#include <boost/dll/import.hpp>
#include <boost/dll/shared_library.hpp>

namespace fast_linalg {

  typedef lapack_int(LAPACKE_dpftrf_t)(int, char, char,
    lapack_int, double *);
  typedef lapack_int(LAPACKE_spftrf_t)(int, char, char,
    lapack_int, float *);
  typedef lapack_int(LAPACKE_dsfrk_t) (int, char, char, char,
    lapack_int, lapack_int, double, const double *,
    lapack_int, double, double *);
  typedef lapack_int(LAPACKE_ssfrk_t) (int, char, char, char,
    lapack_int, lapack_int, double, const float *,
    lapack_int, double, float *);
  typedef lapack_int(LAPACKE_dtfttp_t) (int, char, char,
    lapack_int, const double *, double *);
  typedef lapack_int(LAPACKE_stfttp_t) (int, char, char,
    lapack_int, const float *, float *);
  typedef lapack_int(LAPACKE_dtpttf_t) (int, char, char,
    lapack_int, const double *, double *);
  typedef lapack_int(LAPACKE_stpttf_t) (int, char, char,
    lapack_int, const float *, float *);
  typedef lapack_int(LAPACKE_dpftri_t)(int, char, char,
    lapack_int, double *);
  typedef lapack_int(LAPACKE_spftri_t)(int, char, char,
    lapack_int, float *);
  typedef void (cblas_ssyr_t)(int, int, int, float, const float *,
    int, float *, int);
  typedef void (cblas_dsyr_t)(int, int, int, double, const double *,
    int, double *, int);
  typedef void (cblas_ssyrk_t)(int, int, int,
    int, int, float, float *, int, float, float *, int);
  typedef void (cblas_dsyrk_t)(int, int, int,
    int, int, double, double *, int, double, double *, int);

  typedef int (openblas_get_num_threads_t)();
  typedef int (openblas_get_num_procs_t)();
  typedef void (openblas_set_num_threads_t)(int);
  typedef char* (openblas_get_corename_t)();
  typedef char* (openblas_get_config_t)();

  namespace dll = boost::dll;

  class Wrapper {
    dll::shared_library lib;
    void init_(std::string lib_name) {
      lib = dll::shared_library(lib_name, dll::load_mode::search_system_folders);
      LAPACKE_dpftrf = &lib.get<LAPACKE_dpftrf_t>("LAPACKE_dpftrf");
      LAPACKE_spftrf = &lib.get<LAPACKE_spftrf_t>("LAPACKE_spftrf");
      LAPACKE_dsfrk = &lib.get<LAPACKE_dsfrk_t>("LAPACKE_dsfrk");
      LAPACKE_ssfrk = &lib.get<LAPACKE_ssfrk_t>("LAPACKE_ssfrk");
      LAPACKE_dtfttp = &lib.get<LAPACKE_dtfttp_t>("LAPACKE_dtfttp");
      LAPACKE_stfttp = &lib.get<LAPACKE_stfttp_t>("LAPACKE_stfttp");
      LAPACKE_dtpttf = &lib.get<LAPACKE_dtpttf_t>("LAPACKE_dtpttf");
      LAPACKE_stpttf = &lib.get<LAPACKE_stpttf_t>("LAPACKE_stpttf");
      LAPACKE_dpftri = &lib.get<LAPACKE_dpftri_t>("LAPACKE_dpftri");
      LAPACKE_spftri = &lib.get<LAPACKE_spftri_t>("LAPACKE_spftri");

      openblas_get_num_threads = &lib.get<openblas_get_num_threads_t>(
        "openblas_get_num_threads");
      openblas_get_num_procs = &lib.get<openblas_get_num_procs_t>(
        "openblas_get_num_procs");
      openblas_set_num_threads = &lib.get<openblas_set_num_threads_t>(
        "openblas_set_num_threads");
      openblas_get_corename = &lib.get<openblas_get_corename_t>(
        "openblas_get_corename");
      openblas_get_config = &lib.get<openblas_get_config_t>(
        "openblas_get_config");

      cblas_ssyr = &lib.get<cblas_ssyr_t>("cblas_ssyr");
      cblas_dsyr = &lib.get<cblas_dsyr_t>("cblas_dsyr");
      cblas_ssyrk = &lib.get<cblas_ssyrk_t>("cblas_ssyrk");
      cblas_dsyrk = &lib.get<cblas_dsyrk_t>("cblas_dsyrk");
    }
    void init_() {
      LAPACKE_dpftrf = 0;
      LAPACKE_spftrf = 0;
      LAPACKE_dsfrk = 0;
      LAPACKE_ssfrk = 0;
      LAPACKE_dtfttp = 0;
      LAPACKE_stfttp = 0;
      LAPACKE_dtpttf = 0;
      LAPACKE_stpttf = 0;
      LAPACKE_dpftri = 0;
      LAPACKE_spftri = 0;

      openblas_get_num_threads = 0;
      openblas_get_num_procs = 0;
      openblas_set_num_threads = 0;
      openblas_get_corename = 0;
      openblas_get_config = 0;

      cblas_ssyr = 0;
      cblas_dsyr = 0;
      cblas_ssyrk = 0;
      cblas_dsyrk = 0;
    }
  public:
    LAPACKE_dpftrf_t *LAPACKE_dpftrf;
    LAPACKE_spftrf_t *LAPACKE_spftrf;
    LAPACKE_dsfrk_t *LAPACKE_dsfrk;
    LAPACKE_ssfrk_t *LAPACKE_ssfrk;
    LAPACKE_dtfttp_t *LAPACKE_dtfttp;
    LAPACKE_stfttp_t *LAPACKE_stfttp;
    LAPACKE_dtpttf_t *LAPACKE_dtpttf;
    LAPACKE_stpttf_t *LAPACKE_stpttf;
    LAPACKE_dpftri_t *LAPACKE_dpftri;
    LAPACKE_spftri_t *LAPACKE_spftri;
    openblas_get_num_threads_t *openblas_get_num_threads;
    openblas_get_num_procs_t *openblas_get_num_procs;
    openblas_set_num_threads_t *openblas_set_num_threads;
    openblas_get_corename_t *openblas_get_corename;
    openblas_get_config_t *openblas_get_config;

    cblas_ssyr_t *cblas_ssyr;
    cblas_dsyr_t *cblas_dsyr;
    cblas_ssyrk_t *cblas_ssyrk;
    cblas_dsyrk_t *cblas_dsyrk;

    static Wrapper &instance() {
      SCITBX_ASSERT(initialised());
      return *instance_();
    }

    Wrapper() {
      SCITBX_ASSERT(instance_() == 0);
      init_();
      instance_() = this;
    }

    void initialise(const std::string &lib_name) {
      SCITBX_ASSERT(instance_() != 0);
      init_(lib_name);
    }

    static bool initialised() {
      Wrapper *w = instance_();
      return w != 0 && w->LAPACKE_dpftrf != 0;
    }
    // for internal use only!
    static Wrapper *&instance_() {
      static Wrapper *inst = 0;
      return inst;
    }
  };
}

using namespace fast_linalg;

bool is_fast_linalg_initialised() {
  return Wrapper::initialised();
}
//............................................................................
void initialise_fast_linalg(const std::string &lib_name) {
  if (Wrapper::instance_() == 0) {
    new Wrapper();
  }
  Wrapper::instance_()->initialise(lib_name);
}
//............................................................................
void finalise_fast_linalg() {
  if (Wrapper::instance_() != 0) {
    delete Wrapper::instance_();
    Wrapper::instance_() = 0;
  }
}
//............................................................................
//............................................................................
lapack_int lapack_pftrf(int matrix_order, char transr, char uplo,
    lapack_int n, double* a)
  {
  return (*Wrapper::instance().LAPACKE_dpftrf)(
    matrix_order, transr, uplo, n, a);
}
//............................................................................
lapack_int lapack_spftrf(int matrix_order, char transr, char uplo,
    lapack_int n, float* a)
{
  return (*Wrapper::instance().LAPACKE_spftrf)(
    matrix_order, transr, uplo, n, a);
}
//............................................................................
//............................................................................
lapack_int lapack_sfrk(int matrix_order, char transr, char uplo, char trans,
    lapack_int n, lapack_int k, double alpha,
    const double* a, lapack_int lda, double beta,
    double* c)
{
  return (*Wrapper::instance().LAPACKE_dsfrk)(
    matrix_order, transr, uplo, trans,
    n, k, alpha,
    a, lda, beta,
    c);
}
//............................................................................
lapack_int lapack_ssfrk(int matrix_order, char transr, char uplo, char trans,
    lapack_int n, lapack_int k, float alpha,
    const float* a, lapack_int lda, float beta,
    float* c)
{
  return (*Wrapper::instance().LAPACKE_ssfrk)(
    matrix_order, transr, uplo, trans,
    n, k, alpha,
    a, lda, beta,
    c);
}
//............................................................................
//............................................................................
lapack_int lapack_tfttp(int matrix_order, char transr, char uplo,
    lapack_int n, const double* arf, double* ap)
{
  return (*Wrapper::instance().LAPACKE_dtfttp)(
    matrix_order, transr, uplo, n, arf, ap);
}
//............................................................................
lapack_int lapack_stfttp(int matrix_order, char transr, char uplo,
    lapack_int n, const float* arf, float* ap)
{
  return (*Wrapper::instance().LAPACKE_stfttp)(
    matrix_order, transr, uplo, n, arf, ap);
}
//............................................................................
//............................................................................
lapack_int lapack_tpttf(int matrix_order, char transr, char uplo,
    lapack_int n, const double* ap, double* arf)
{
  return (*Wrapper::instance().LAPACKE_dtpttf)(
    matrix_order, transr, uplo, n, ap, arf);
}
//............................................................................
lapack_int lapack_stpttf(int matrix_order, char transr, char uplo,
    lapack_int n, const float* ap, float* arf)
{
  return (*Wrapper::instance().LAPACKE_stpttf)(
    matrix_order, transr, uplo, n, ap, arf);
}
//............................................................................
//............................................................................
lapack_int lapack_pftri(int matrix_order, char transr, char uplo,
    lapack_int n, double* a)
{
  return (*Wrapper::instance().LAPACKE_dpftri)(
    matrix_order, transr, uplo, n, a);
}
//............................................................................
lapack_int lapack_spftri(int matrix_order, char transr, char uplo,
    lapack_int n, float* a)
{
  return (*Wrapper::instance().LAPACKE_spftri)(
    matrix_order, transr, uplo, n, a);
}
//............................................................................
//............................................................................
void cblas_ssyr(int Order, int Uplo, int N, float alpha,
  const float *X, int incX, float *A, int lda)
{
  return (*Wrapper::instance().cblas_ssyr)(Order, Uplo, N, alpha,
    X, incX, A, lda);
}
//............................................................................
void cblas_dsyr(int Order, int Uplo, int N, double alpha,
  const double *X, int incX, double *A, int lda)
{
  return (*Wrapper::instance().cblas_dsyr)(Order, Uplo, N, alpha,
    X, incX, A, lda);
}
//............................................................................
void cblas_ssyrk(int Order, int Uplo, int Trans, int N, int K,
  float alpha, float *A, int lda, float beta, float *C, int ldc)
{
  return (*Wrapper::instance().cblas_ssyrk)(Order, Uplo, Trans, N, K, alpha,
    A, lda, beta, C, ldc);
}
//............................................................................
void cblas_dsyrk(int Order, int Uplo, int Trans, int N, int K,
  double alpha, double *A, int lda, double beta, double *C, int ldc)
{
  return (*Wrapper::instance().cblas_dsyrk)(Order, Uplo, Trans, N, K, alpha,
    A, lda, beta, C, ldc);
}
//............................................................................
//............................................................................
int openblas_get_num_threads() {
  return (*Wrapper::instance().openblas_get_num_threads)();
}
//............................................................................
int openblas_get_num_procs() {
  return (*Wrapper::instance().openblas_get_num_procs)();
}
//............................................................................
void openblas_set_num_threads(int n) {
  (*Wrapper::instance().openblas_set_num_threads)(n);
}
//............................................................................
char* openblas_get_corename() {
  return (*Wrapper::instance().openblas_get_corename)();
}
//............................................................................
char* openblas_get_config() {
  return (*Wrapper::instance().openblas_get_config)();
}
//............................................................................
