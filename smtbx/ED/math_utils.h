#pragma once
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/versa_matrix.h>
#include <smtbx/error.h>
#include <smtbx/import_scitbx_af.h>
#include <smtbx/ED/ed.h>

namespace smtbx {  namespace ED
{
  using namespace cctbx;

  template <typename FloatType>
  struct math_utils {
    ED_UTIL_TYPEDEFS;

    static void multiply_diagonal_inplace(
      af::shared<complex_t>& d,
      af::shared<complex_t> const& c)
    {
      for (size_t i = 0; i < d.size(); i++) {
        d[i] *= c[i];
      }
    }

    static af::shared<FloatType> multiply_diagonal(
      af::shared<complex_t> const& d,
      af::shared<complex_t> const& c)
    {
      af::shared<FloatType> rv(d.size());
      for (size_t i = 0; i < d.size(); i++) {
        rv[i] = d[i] * c[i];
      }
      return rv;
    }

    static cmat_t multiply_diagonal(
      cmat_t const& m,
      af::shared<complex_t> const& d)
    {
      size_t sz = m.accessor().n_columns();
      cmat_t rv(af::mat_grid(sz, sz));
      for (size_t i = 0; i < sz; i++) {
        for (size_t j = 0; j < sz; j++) {
          rv(i, j) = m(i, j) * d[j];
        }
      }
      return rv;
    }

    static void add_mat(cmat_t& A, const cmat_t& B, FloatType sig = 1) {
      size_t sz = A.accessor().n_columns();
      for (size_t i = 0; i < sz; i++) {
        for (size_t j = 0; j < sz; j++) {
          A(i, j) += B(i, j) * sig;
        }
      }
    }

    static cmat_t multiply_diagonal(
      af::shared<complex_t> const& d,
      cmat_t const& m)
    {
      size_t sz = m.accessor().n_columns();
      cmat_t rv(af::mat_grid(sz, sz));
      for (size_t i = 0; i < sz; i++) {
        for (size_t j = 0; j < sz; j++) {
          rv(i, j) = m(i, j) * d[i];
        }
      }
      return rv;
    }

    // derivatites of A = V * D * V^-1
    template <typename ev_t>
    static cmat_t calc_dS_dx(
      // derivatives of the Eigen matrix by params
      cmat_t const& S,
      // diagonal matrix
      af::shared<complex_t> const& D,
      // eigen values
      af::shared<ev_t> const& L,
      // eigen vectors
      cmat_t const& V,
      // inverse of eigen vectors
      cmat_t const& V_1,
      // derivatives of the diagonal matrix by eigen values
      af::shared<complex_t> const& dD_dL)
    {
      size_t sz = S.accessor().n_columns();
      cmat_t X = af::matrix_multiply(V_1.const_ref(), S.const_ref());
      cmat_t dV = af::matrix_multiply(X.const_ref(), V.const_ref());
      af::shared<complex_t> dD(sz);
      for (size_t i = 0; i < sz; i++) {
        dD[i] = dV(i, i) * dD_dL[i];
      }

      for (size_t i = 0; i < sz; i++) {
        dV(i, i) = 0;
        for (size_t j = i + 1; j < sz; j++) {
          complex_t k = 1. / (L[j] - L[i]);
          dV(i, j) *= k;
          dV(j, i) *= -k;
        }
      }
      dV = af::matrix_multiply(V.const_ref(), dV.const_ref());
      X = multiply_diagonal(D, V_1);
      cmat_t res = af::matrix_multiply(dV.const_ref(), X.const_ref());

      add_mat(res,
        af::matrix_multiply(
          multiply_diagonal(V, dD).const_ref(), V_1.const_ref())
      );

      add_mat(res,
        af::matrix_multiply(
          af::matrix_multiply(
            af::matrix_multiply(V.const_ref(), X.const_ref()).const_ref(),
            dV.const_ref()).const_ref(),
          V_1.const_ref()),
        -1
      );
      return res;
    }

    /* derivatites of A = V * D * V^-1 *col(1,0,0,..)
    tested agains the first col of the full procedure above
    */
    template <typename ev_t>
    static af::shared<complex_t> calc_dS_dx_1(
      // derivatives of the Eigen matrix by params
      cmat_t const& S,
      // diagonal matrix
      af::shared<complex_t> const& D,
      // eigen values
      af::shared<ev_t> const& L,
      // eigen vectors
      cmat_t const& V,
      // inverse of eigen vectors
      cmat_t const& V_1,
      // derivatives of the diagonal matrix by eigen values
      af::shared<complex_t> const& dD_dL)
    {
      size_t sz = S.accessor().n_columns();
      cmat_t X = af::matrix_multiply(V_1.const_ref(), S.const_ref());
      cmat_t dV = af::matrix_multiply(X.const_ref(), V.const_ref());

      af::shared<complex_t> dD(sz);
      for (size_t i = 0; i < sz; i++) {
        dD[i] = dV(i, i) * dD_dL[i];
      }

      for (size_t i = 0; i < sz; i++) {
        dV(i, i) = 0;
        for (size_t j = i + 1; j < sz; j++) {
          complex_t k = 1. / (L[j] - L[i]);
          dV(i, j) *= k;
          dV(j, i) *= -k;
        }
      }
      dV = af::matrix_multiply(V.const_ref(), dV.const_ref());

      // first col of inverse of eigen vectors
      af::shared<complex_t> V_1_1(sz), tmp(sz);
      for (size_t i = 0; i < sz; i++) {
        V_1_1[i] = V_1(i, 0);
      }

      // evaluate first term dV*D*V_1_1
      for (size_t i = 0; i < sz; i++) {
        tmp[i] = D[i] * V_1_1[i];
      }
      af::shared<complex_t> d_V = af::matrix_multiply(dV.const_ref(), tmp.const_ref());

      // evaluate second term V*d_D*V_1_1
      for (size_t i = 0; i < sz; i++) {
        tmp[i] = dD[i] * V_1_1[i];
      }
      af::shared<complex_t> d_L = af::matrix_multiply(V.const_ref(), tmp.const_ref());

      // evaluate third term V*D*(V_1*dV*V_1_1)
      tmp = af::matrix_multiply(dV.const_ref(), V_1_1.const_ref());
      tmp = af::matrix_multiply(V_1.const_ref(), tmp.const_ref());
      for (size_t i = 0; i < sz; i++) {
        tmp[i] *= D[i];
      }
      af::shared<complex_t> d_V_1 = af::matrix_multiply(V.const_ref(), tmp.const_ref());
      // sum up
      for (size_t i = 0; i < sz; i++) {
        d_V[i] = d_V[i] + d_L[i] - d_V_1[i];
      }
      return d_V;
    }

    // Assumes A(0,0)=0, replaces A with column eigen vecs
    //https://quantumcomputing.stackexchange.com/questions/22222/how-to-find-the-eigenstates-of-a-general-2-times-2-hermitian-matrix
    static void two_beam_eigen(complex_t* A,
      FloatType* ev)
    {
      FloatType h11 = A[3].real() / 2;
      FloatType A01_sq = std::norm(A[1]);
      FloatType s = std::sqrt(h11 * h11 + A01_sq);
      ev[0] = h11 + s;
      ev[1] = h11 - s;
      FloatType v1l = std::sqrt(A01_sq + ev[0] * ev[0]);
      FloatType v2l = std::sqrt(A01_sq + ev[1] * ev[1]);
      complex_t A01 = A[1];
      A[0] = A01 / v1l;  A[1] = ev[0] / v1l;
      A[2] = A01 / v2l;  A[3] = ev[1] / v2l;
    }


  }; //struct smtbx::ED::math_utils
}} // namespace smtbx::ED