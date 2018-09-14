#pragma once

namespace fast_linalg {
  const int CblasRowMajor = 101,
    CblasColMajor = 102;
  const int CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113,
    CblasConjNoTrans = 114;
  const int CblasUpper = 121,
    CblasLower = 122;
  const int CblasNonUnit = 131,
    CblasUnit = 132;
  const int CblasLeft = 141,
    CblasRight = 142;
}
#if defined(USE_FAST_LINALG)
#include <boost/config.hpp>

#ifndef fast_linalg_api
#define fast_linalg_api BOOST_SYMBOL_EXPORT
#endif

extern "C" {
  fast_linalg_api int openblas_get_num_threads();
  fast_linalg_api int openblas_get_num_procs();
  fast_linalg_api void openblas_set_num_threads(int);
  fast_linalg_api char* openblas_get_corename();
  fast_linalg_api char* openblas_get_config();
};
#endif // USE_FAST_LINALG
