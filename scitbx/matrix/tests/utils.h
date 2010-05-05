#ifndef SCITBX_MATRIX_TESTS_UTILS_H
#define SCITBX_MATRIX_TESTS_UTILS_H

#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/accessors/packed_matrix.h>

namespace scitbx { namespace matrix {

typedef af::c_grid<2> dim;
typedef af::versa<double, dim> matrix_t;
typedef af::ref<double, dim> matrix_ref_t;
typedef af::const_ref<double, dim> matrix_const_ref_t;
typedef af::versa<double, af::packed_u_accessor> symmetric_matrix_packed_u_t;
typedef af::ref<double, af::packed_u_accessor> symmetric_matrix_packed_u_ref_t;
typedef af::versa<double, af::packed_l_accessor> symmetric_matrix_packed_l_t;
typedef af::ref<double, af::packed_l_accessor> symmetric_matrix_packed_l_ref_t;
typedef af::shared<double> vec_t;
typedef af::ref<double> vec_ref_t;

matrix_t product_U_M_VT(matrix_const_ref_t const &u,
                        matrix_const_ref_t const &m,
                        matrix_const_ref_t const &v)
{
  matrix_t u_m = af::matrix_multiply(u, m);
  matrix_t vt = af::matrix_transpose(v);
  return af::matrix_multiply(u_m.const_ref(), vt.const_ref());
}

matrix_t product_UT_M_V(matrix_const_ref_t const &u,
                        matrix_const_ref_t const &m,
                        matrix_const_ref_t const &v)
{
  matrix_t ut = af::matrix_transpose(u);
  matrix_t m_v = af::matrix_multiply(m, v);
  return af::matrix_multiply(ut.const_ref(), m_v.const_ref());
}

}} // scitbx::matrix

#endif // GUARD
