#pragma once
#include <cctbx/miller.h>
#include <cctbx/miller/index_generator.h>
#include <cctbx/miller/lookup_utils.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/versa_matrix.h>
#include <smtbx/error.h>
#include <smtbx/import_scitbx_af.h>
#include <fast_linalg/lapacke.h>

namespace smtbx { namespace ED {
  using namespace cctbx;

  template <typename FloatType>
  struct utils {
    typedef std::complex<FloatType> complex_t;
    typedef scitbx::vec3<FloatType> cart_t;
    typedef scitbx::mat3<FloatType> mat3_t;
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;

    static void build_Ug_matrix(
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<complex_t> const& Fcs_k,
      lookup_t const& mi_lookup,
      af::shared<miller::index<> > indices, // implicit {0,0,0} at 0
      FloatType Fc2Ug)
    {
      const size_t n_beams = indices.size() + 1; // g0+
      A.resize(af::mat_grid(n_beams, n_beams));
      A(0, 0) = 0;
      for (size_t i = 1; i < n_beams; i++) {
        miller::index<> h_i = indices[i - 1];
        int ii = mi_lookup.find_hkl(h_i);
        complex_t Fc_i = ii != -1 ? Fcs_k[ii] : 0;
        // h_i - (0,0,0)
        A(i, 0) = Fc2Ug * Fc_i;
        // (0,0,0) - h_i
        A(0, i) = std::conj(A(i, 0));
        A(i, i) = 0;
        for (size_t j = i + 1; j < n_beams; j++) {
          miller::index<> h_j = indices[j - 1];
          int i_m_j = mi_lookup.find_hkl(h_i - h_j);
          complex_t Fc = 0;
          if (i_m_j < 0) {
            int j_m_i = mi_lookup.find_hkl(h_j - h_i);
            if (j_m_i >= 0) {
              Fc = std::conj(Fcs_k[j_m_i]);
            }
          }
          else {
            Fc = Fcs_k[i_m_j];
          }
          A(i, j) = Fc2Ug * Fc;
          A(j, i) = std::conj(A(i, j));
        }
      }
    }

    // considers the original matrix Hermitian
    static void modify_Ug_matrix(
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<miller::index<> > const& s_indices, // matrix beams
      af::shared<complex_t> const& Fcs_k,
      lookup_t const& mi_lookup,
      af::shared<miller::index<> > const& w_indices, // weak beams for pertubation
      af::shared<FloatType> const& Exitations, // 2KSg for w_indices
      FloatType Fc2Ug)
    {
      const size_t n_beams = A.accessor().n_columns();
      FloatType k = Fc2Ug * Fc2Ug;
      for (size_t i = 1; i < n_beams; i++) {
        miller::index<> h_i = s_indices[i - 1];
        complex_t d_mod = 0, u_mod = 0;
        for (size_t j = 0; j < w_indices.size(); j++) {
          miller::index<> h_j = w_indices[j];
          int i_j = mi_lookup.find_hkl(h_j);
          complex_t Fc_j = i_j != -1 ? Fcs_k[i_j] : 0;
          int i_m_j = mi_lookup.find_hkl(h_i - h_j);
          complex_t Fc_i_m_j = 0;
          if (i_m_j < 0) {
            int j_m_i = mi_lookup.find_hkl(h_j - h_i);
            if (j_m_i >= 0) {
              Fc_i_m_j = std::conj(Fcs_k[j_m_i]);
            }
          }
          else {
            Fc_i_m_j = Fcs_k[i_m_j];
          }
          u_mod += Fc_i_m_j * Fc_j / Exitations[j];
          d_mod += Fc_i_m_j * std::conj(Fc_i_m_j) / Exitations[j];
        }
        A(i, i) -= k * d_mod;
        A(i, 0) -= k * u_mod;
      }
    }

    /* Acta Cryst. (2013). A69, 171–188 */
    static void build_eigen_matrix_2013(
      af::versa<complex_t, af::mat_grid>& A, // Ug matrix
      af::shared<miller::index<> > indices, // implicit {0,0,0} at 0
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N,
      af::shared<FloatType>& ExpDen,
      FloatType Fc2Ug)
    {
      const FloatType Kn = N * K, Kl = K.length();
      const size_t n_beams = indices.size() + 1; // g0+
      ExpDen.resize(n_beams);
      ExpDen[0] = Kn; //for g0
      for (size_t i = 1; i < n_beams; i++) {
        miller::index<> h = indices[i - 1];
        cart_t K_g = K + RMf * cart_t(h[0], h[1], h[2]);
        FloatType s_2k = Kl * Kl - K_g.length_sq();
        A(i, i) += s_2k;
        ExpDen[i] = K_g * N;
      }
    }

    /* Acta Cryst. (2015). A71, 235–244 */
    static void build_eigen_matrix_2015(
      af::versa<complex_t, af::mat_grid>& A, // Ug matrix
      af::shared<miller::index<> > indices, // implicit {0,0,0} at 0
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N,
      af::shared<FloatType>& M,
      FloatType Fc2Ug)
    {
      const FloatType Kn = N * K,
        Kl = K.length();
      const size_t n_beams = indices.size() + 1; // g0+
      M.resize(n_beams);
      M[0] = 1; //for g0
      af::shared<cart_t> gs(n_beams);
      af::shared<FloatType> dens(n_beams);
      dens[0] = 1;
      for (size_t i = 1; i < n_beams; i++) {
        miller::index<> h = indices[i - 1];
        cart_t g = RMf * cart_t(h[0], h[1], h[2]);
        gs[i] = g;
        dens[i] = std::sqrt(1. / (1 + g * N / Kn));
      }
      for (size_t i = 1; i < n_beams; i++) {
        FloatType s_2k = (Kl * Kl - (K + gs[i]).length_sq());
        FloatType i_den = dens[i];
        A(0, i) *= i_den;
        A(i, i) *= i_den;
        A(i, i) += s_2k * i_den * i_den;
        M[i] = i_den;
        for (size_t j = i + 1; j < n_beams; j++) {
          FloatType den = dens[j] * i_den;
          A(i, j) *= den;
          A(j, i) *= den;
        }
      }
    }

    /* J. Appl. Cryst. (2022). 55 */
    static void build_eigen_matrix_recipro(
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<miller::index<> > const &indices,
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N,
      af::shared<FloatType>& Pgs,
      FloatType Fc2Ug)
    {
      // Projection of K onto normal of frame normal, K*frame.normal
      const FloatType Kn = N * K;
      const size_t n_beams = indices.size() + 1;
      A.resize(af::mat_grid(n_beams, n_beams));
      Pgs.resize(n_beams);
      for (size_t i = 0; i < n_beams; i++) {
        miller::index<> h_i = i == 0 ? miller::index<>(0,0,0) : indices[i-1];
        cart_t g_i = RMf * cart_t(h_i[0], h_i[1], h_i[2]);
        FloatType Pg = 2 * (N * (K + g_i));
        Pgs[i] = Pg;
        for (size_t j = 0; j < n_beams; j++) {
          if (i == j) {
            A(i, i) += -(g_i * (K + K + g_i));
          }
          A(i, j) /= Pg;
        }
      }
    }

    static void build_eigen_matrix_modified(
      af::versa<complex_t, af::mat_grid>& A, // Ug matrix modified
      cart_t const& K,
      cart_t const& N,
      af::shared<cart_t> const& gs,
      af::shared<FloatType> const& excitation_errors,
      af::shared<FloatType>& ExpDen)
    {
      const FloatType Kn = N * K, Kl = K.length();
      const size_t n_beams = A.accessor().n_columns();
      ExpDen.resize(n_beams);
      ExpDen[0] = Kn; //for g0
      for (size_t i = 1; i < n_beams; i++) {
        A(i, i) += excitation_errors[i - 1];
        ExpDen[i] = (K + gs[i - 1]) * N;
      }
    }

    static af::shared<complex_t> calc_amps_2013(
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<FloatType> const& ExpDen,
      FloatType thickness,
      size_t num)
    {
      using namespace fast_linalg;
      const size_t n_beams = A.accessor().n_columns();
      af::shared<FloatType> ev(n_beams);
      // heev replaces A with column-wise eigenvectors
      lapack_int info = heev(LAPACK_ROW_MAJOR, 'V', LAPACK_UPPER, n_beams,
        A.begin(), n_beams, ev.begin());
      SMTBX_ASSERT(!info)(info);
      const complex_t exp_k(0, scitbx::constants::pi * thickness);
      af::shared<complex_t> im(n_beams);
      for (size_t i = 0; i < n_beams; i++) {
        im[i] = std::exp(ev[i] * exp_k / ExpDen[i]) * std::conj(A(0, i));
      }
      af::shared<complex_t> res = af::matrix_multiply(A.const_ref(), im.const_ref());
      af::shared<complex_t> rv;
      for (size_t i = 1; i < n_beams; i++) {
        rv.push_back(res[i]);
        if (rv.size() >= num) {
          break;
        }
      }
      return rv;
    }

    static af::shared<complex_t> calc_amps_2015(
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<FloatType> const& M,
      FloatType thickness,
      FloatType Kn,
      size_t num)
    {
      using namespace fast_linalg;
      const size_t n_beams = A.accessor().n_columns();
      af::shared<FloatType> ev(n_beams);
      // heev replaces A with column-wise eigenvectors
      lapack_int info = heev(LAPACK_ROW_MAJOR, 'V', LAPACK_UPPER, n_beams,
        A.begin(), n_beams, ev.begin());
      SMTBX_ASSERT(!info)(info);
      af::shared<complex_t> im(n_beams);
      const complex_t exp_k(0, scitbx::constants::pi * thickness / Kn);
      for (size_t i = 0; i < n_beams; i++) {
        im[i] = std::exp(ev[i] * exp_k) * std::conj(A(0, i));
      }
      af::shared<complex_t> res = af::matrix_multiply(A.const_ref(), im.const_ref());
      af::shared<complex_t> rv;
      for (size_t i = 1; i < n_beams; i++) {
        rv.push_back(res[i] * M[i]);
        if (rv.size() >= num) {
          break;
        }
      }
      return rv;
    }

    static af::shared<complex_t> calc_amps_recipro(
      af::versa<complex_t, af::mat_grid>& A,
      const af::shared<FloatType>& Pgs,
      FloatType Kvac,
      FloatType thickness,
      size_t num)
    {
      using namespace fast_linalg;
      const size_t n_beams = A.accessor().n_columns();

      // eigenvalues, this is for generic matrix
      af::shared<complex_t> ev(n_beams);
      // right eigenvectors
      af::versa<complex_t, af::c_grid<2> > B(af::c_grid<2>(n_beams, n_beams));
      lapack_int info = geev(LAPACK_ROW_MAJOR, 'N', 'V', n_beams,
        A.begin(), n_beams, ev.begin(), 0, n_beams, B.begin(), n_beams);
      SMTBX_ASSERT(!info)(info);
      af::versa<complex_t, af::mat_grid> Bi = B.deep_copy();
      // invert Bi now
      {
        af::shared<lapack_int> pivots(n_beams);
        info = getrf(LAPACK_ROW_MAJOR, n_beams, n_beams, Bi.begin(),
          n_beams, pivots.begin());
        SMTBX_ASSERT(!info)(info);
        info = getri(LAPACK_ROW_MAJOR, n_beams, Bi.begin(),
          n_beams, pivots.begin());
        SMTBX_ASSERT(!info)(info);
      }
      af::shared<complex_t> im(n_beams);
      const complex_t exp_k(0, 2*scitbx::constants::pi * thickness);
      const complex_t exp_k1(0, scitbx::constants::pi * thickness);
      // apply diagonal matrix on the left
      for (size_t i = 0; i < n_beams; i++) {
        complex_t p = std::exp(exp_k1 * (Pgs[i] + 2 * Kvac));
        for (size_t j = 0; j < n_beams; j++) {
          B(i, j) *= p;
        }
      }
      for (size_t i = 0; i < n_beams; i++) {
        im[i] = std::exp(ev[i] * exp_k) * Bi(i, 0);
      }
      af::shared<complex_t> res = af::matrix_multiply(B.const_ref(), im.const_ref());
      af::shared<complex_t> rv;
      for (size_t i = 1; i < n_beams; i++) {
        rv.push_back(res[i]);
        if (rv.size() >= num) {
          break;
        }
      }
      return rv;
    }

    static af::shared<complex_t> calc_amps_modified(
      af::versa<complex_t, af::mat_grid>& A,
      const af::shared<FloatType>& ExpDen,
      FloatType thickness,
      size_t num)
    {
      using namespace fast_linalg;
      const size_t n_beams = A.accessor().n_columns();

      // eigenvalues, this is for generic matrix
      af::shared<complex_t> ev(n_beams);
      // right eigenvectors
      af::versa<complex_t, af::c_grid<2> > B(af::c_grid<2>(n_beams, n_beams));
      lapack_int info = geev(LAPACK_ROW_MAJOR, 'N', 'V', n_beams,
        A.begin(), n_beams, ev.begin(), 0, n_beams, B.begin(), n_beams);
      SMTBX_ASSERT(!info)(info);
      af::versa<complex_t, af::mat_grid> Bi = B.deep_copy();
      // invert Bi now
      {
        af::shared<lapack_int> pivots(n_beams);
        info = getrf(LAPACK_ROW_MAJOR, n_beams, n_beams, Bi.begin(),
          n_beams, pivots.begin());
        SMTBX_ASSERT(!info)(info);
        info = getri(LAPACK_ROW_MAJOR, n_beams, Bi.begin(),
          n_beams, pivots.begin());
        SMTBX_ASSERT(!info)(info);
      }

      const complex_t exp_k(0, scitbx::constants::pi * thickness);
      af::shared<complex_t> im(n_beams);
      for (size_t i = 0; i < n_beams; i++) {
        im[i] = std::exp(ev[i] * exp_k / ExpDen[i]) * Bi(i, 0);
      }
      af::shared<complex_t> res = af::matrix_multiply(B.const_ref(), im.const_ref());
      af::shared<complex_t> rv;
      for (size_t i = 1; i < n_beams; i++) {
        rv.push_back(res[i]);
        if (rv.size() >= num) {
          break;
        }
      }
      return rv;
    }

    // Assumes A(0,0)=0, replaces A with column egein vecs
    //https://quantumcomputing.stackexchange.com/questions/22222/how-to-find-the-eigenstates-of-a-general-2-times-2-hermitian-matrix
    static void two_beam_eigen(af::versa<complex_t, af::mat_grid> &A,
      af::shared<FloatType> &ev)
    {
      FloatType h11 = A(1, 1).real() / 2;
      FloatType s = std::sqrt(h11 * h11 + std::norm(A(0, 1)));
      ev[0] = h11 + s;
      ev[1] = h11 - s;
      FloatType v1l = std::sqrt(2 * s * ev[0]);
      FloatType v2l = std::sqrt(-2 * s * ev[1]);
      complex_t A01 = A(0, 1);
      A(0, 0) = A01 / v1l;  A(0, 1) = ev[0] / v1l;
      A(1, 0) = A01 / v2l;  A(1, 1) = ev[1] / v2l;
    }

    static complex_t calc_amp_2beam(
      const miller::index<> &h, const complex_t Ug,
      FloatType thickness,
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N)
    {
      using namespace fast_linalg;
      const FloatType Kn = N * K, Kl = K.length();
      cart_t K_g = K + RMf * cart_t(h[0], h[1], h[2]);
      FloatType s_2k = Kl * Kl - K_g.length_sq();

      af::versa<complex_t, af::mat_grid> A(af::mat_grid(2,2));
      A(0, 0) = 0;
      A(1, 0) = Ug;
      A(0, 1) = std::conj(Ug);
      A(1, 1) = s_2k;
      af::shared<FloatType> ev(2);
      //two_beam_eigen(A, ev);
      // heev replaces A with column-wise eigenvectors
      lapack_int info = heev(LAPACK_ROW_MAJOR, 'V', LAPACK_UPPER, 2,
        A.begin(), 2, ev.begin());
      SMTBX_ASSERT(!info)(info);
      const complex_t exp_k(0, scitbx::constants::pi * thickness);
      af::shared<complex_t> im(2);
      FloatType ExpDen = K_g * N;
      im[0] = std::exp(ev[0] * exp_k / Kn) * std::conj(A(0, 0));
      im[1] = std::exp(ev[1] * exp_k / (K_g * N)) * std::conj(A(0, 1));
      af::shared<complex_t> res = af::matrix_multiply(A.const_ref(), im.const_ref());
      return res[1];
    }

    static bool is_excited_g(cart_t const& g,
      cart_t const& K, FloatType MaxSg, FloatType MaxG, FloatType precession_angle)
    {
      FloatType gl = g.length(), Kl = K.length();
      cart_t Kg = K + g;
      FloatType Sg = std::abs((Kl * Kl - Kg.length_sq()) / (2 * Kl));
      return Sg < MaxSg && gl < MaxG && Sg / (gl * precession_angle) < 1.0;
    }

    static bool is_excited_h(miller::index<> const& index,
      mat3_t const& RMf, // matrix to orthogonalise and rotate into the frame basis
      cart_t const& K, FloatType MaxSg, FloatType MaxG, FloatType precession_angle)
    {
      return is_excited_g(RMf * cart_t(index[0], index[1], index[2]), K,
        MaxSg, MaxG, precession_angle);
    }

    struct ExcitedBeam {
      miller::index<> h;
      cart_t g;
      FloatType w, Sg;
      ExcitedBeam()
        : w(0), Sg(0)
      {}
      ExcitedBeam(miller::index<> const& index, cart_t const& g, FloatType weight,
        FloatType Sg)
        : h(index), g(g),
        w(weight), Sg(Sg)
      {}
      static bool compare(ExcitedBeam const& a, ExcitedBeam const& b) {
        return a.w < b.w;
      }
    };

    /* Generates a list of miller indices for given resolution and space group
    * The indices are will fulfil the MaxSg and MaxG parameters and sorted by
    * FoM from ReciPro
    */
    static af::shared<ExcitedBeam> generate_index_set(
      mat3_t const& RMf, // matrix to orthogonalise and rotate into the frame basis
      cart_t K,
      FloatType min_d,
      FloatType MaxG, FloatType MaxSg,
      uctbx::unit_cell const& unit_cell)
    {
      using namespace cctbx::miller;

      index_generator h_generator(unit_cell, sgtbx::space_group_type("P1"),
        true,
        min_d, true);

      miller::index<> h;
      af::shared<ExcitedBeam> all;
      
      const FloatType max_f_sq = MaxG * MaxG,
        Kl = K.length(),
        Kl_sq = Kl * Kl;
      while (!(h = h_generator.next()).is_zero()) {
        cart_t g = RMf * cart_t(h[0], h[1], h[2]);
        FloatType g_sq = g.length_sq();
        if (g_sq > max_f_sq) {
          continue;
        }
        FloatType Kg_sq = (K + g).length_sq();
        FloatType s = Kl_sq - Kg_sq;
        FloatType Sg = std::abs(s / (2 * Kl));
        if (MaxSg > 0 && Sg > MaxSg) {
          continue;
        }
        FloatType w = s * s * Kg_sq;
        all.push_back(ExcitedBeam(h, g, w, Sg));
      }
      std::sort(all.begin(), all.end(), &ExcitedBeam::compare);
      return all;
    }

  }; //struct smtbx::ED::utils
}} // namespace smtbx::ED