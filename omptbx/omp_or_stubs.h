#ifndef OMPTBX_WRAPPER_H
#define OMPTBX_WRAPPER_H

#if defined(_OPENMP)
# define OMPTBX_HAVE_OMP_H
# include <omp.h>
#else
# define OMPTBX_HAVE_STUBS_H
# include <omptbx/stubs.h>
#endif

#endif // OMPTBX_WRAPPER_H
