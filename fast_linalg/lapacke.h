#pragma once

#if defined(CCTBX_FAST_LINALG_USES_OPENBLAS)
#define LAPACK_COMPLEX_CUSTOM
#include <complex>
#define lapack_complex_float  std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>

#else
#error No LAPACKE has been configured
#endif
