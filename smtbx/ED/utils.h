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

    static bool is_excited_g(cart_t const& g_,
      FloatType Kl, FloatType MaxSg, FloatType MaxG, FloatType precession_angle)
    {
      FloatType gl = g_.length();
      cart_t g(g_[0], g_[1], g_[2] - Kl);
      FloatType Sg = std::abs((Kl * Kl - g.length_sq()) / (2 * Kl));
      return Sg < MaxSg&& gl < MaxG && Sg / (gl * precession_angle) < 1.0;
    }

    static bool is_excited_h(miller::index<> const& index,
      mat3_t const& RMf, // matrix to orthogonalise and rotate into the frame basis
      FloatType Kl, FloatType MaxSg, FloatType MaxG, FloatType precession_angle)
    {
      return is_excited_g(RMf * cart_t(index[0], index[1], index[2]), Kl,
        MaxSg, MaxG, precession_angle);
    }

    struct ExcitedBeam {
      miller::index<> h;
      cart_t g;
      FloatType w, s;
      ExcitedBeam(miller::index<> const& index, cart_t const& g, FloatType weight, FloatType s)
        : h(index), g(g),
        w(weight), s(s)
      {}
      static bool compare(ExcitedBeam const& a, ExcitedBeam const& b) {
        return a.w < b.w;
      }
    };

    /* Generates a list of miller indices for given resolution and space group
    * The indices are sorted by excitation FoM asc
    */
    static af::shared<ExcitedBeam> generate_index_set(
      mat3_t const& RMf, // matrix to orthogonalise and rotate into the frame basis
      FloatType Kl,
      FloatType min_d,
      uctbx::unit_cell const& unit_cell,
      sgtbx::space_group const& space_group, bool anomalous)
    {
      using namespace cctbx::miller;

      index_generator h_generator(unit_cell, space_group.type(), anomalous, min_d);

      miller::index<> h;
      af::shared<ExcitedBeam> all;
      while (!(h = h_generator.next()).is_zero()) {
        cart_t g = RMf * cart_t(h[0], h[1], h[2]);
        g[2] += Kl;
        FloatType s = std::abs(Kl * Kl - g.length_sq());
        FloatType w = s * s * g.length_sq();
        all.push_back(ExcitedBeam(h, g, w, s));
      }
      std::sort(all.begin(), all.end(), &ExcitedBeam::compare);
      return all;
    }

    static af::shared<miller::index<> > filter_index_set(
      const af::shared<ExcitedBeam>& beams,
      FloatType Kl, FloatType MaxSg, FloatType MaxG, FloatType precession_angle)
    {
      af::shared<miller::index<> > rv;
      for (size_t i = 0; i < beams.size(); i++) {
        if (is_excited_g(beams[i].g, Kl, MaxSg, MaxG, precession_angle)) {
          rv.push_back(beams[i].h);
        }
      }
      return rv;
    }

    static af::shared<ExcitedBeam> update_index_set(
      af::shared<ExcitedBeam> &beams,
      mat3_t const& RMf, // matrix to orthogonalise and rotate into the frame basis
      FloatType Kl)
    {
      using namespace cctbx::miller;
      for (size_t i = 0; i < beams.size(); i++) {
        cart_t g = RMf * cart_t(beams[i].h[0], beams[i].h[1], beams[i].h[2]);
        g[2] += Kl;
        FloatType s = std::abs(Kl * Kl - g.length_sq());
        beams[i].w = s * s * g.length_sq();
        beams[i].s = s;
        beams[i].g = g;
      }
      std::sort(beams.begin(), beams.end(), &ExcitedBeam::compare);
      return beams;
    }
  }; //struct smtbx::ED::utils
}} // namespace smtbx::ED