#pragma once
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/versa_matrix.h>
#include <smtbx/error.h>
#include <smtbx/import_scitbx_af.h>
#include<cctbx/miller.h>
#include <cctbx/miller/lookup_utils.h>
#include <scitbx/array_family/selections.h>
#include <fast_linalg/lapacke.h>

#define ED_UTIL_TYPEDEFS                                           \
  typedef std::complex<FloatType> complex_t;                       \
  typedef scitbx::vec3<FloatType> cart_t;                          \
  typedef scitbx::mat3<FloatType> mat3_t;                          \
  typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t; \
  typedef typename af::versa<FloatType, af::mat_grid> mat_t;       \
  typedef typename af::versa<complex_t, af::mat_grid> cmat_t;

namespace smtbx { namespace ED{
  enum {
    /* As in Acta Cryst. (2013). A69, 171–188 */
    DYN_MATRIX_2013 = 0,
    /* As in Acta Cryst. (2015). A71, 235–244 */
    DYN_MATRIX_2015 = 1,
    /* As in Electron Microscopy of Thin Crystals by Hirsch, Peter B., etc., et al
    (ISBN: 9780882753768)
    */
    DYN_MATRIX_DEFAULT = 10
  };
}}
