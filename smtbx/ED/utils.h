#pragma once
#include <cctbx/miller.h>
#include <cctbx/miller/index_generator.h>
#include <smtbx/ED/math_utils.h>

namespace smtbx { namespace ED
{
  using namespace cctbx;

  template <typename FloatType>
  struct utils {
    ED_UTIL_TYPEDEFS;

    static void build_Ug_matrix(
      cmat_t& A,
      const af::shared<complex_t> &Fcs_k,
      const lookup_t &mi_lookup,
      const af::shared<miller::index<> > &indices) // implicit {0,0,0} at 0
    {
      const size_t n_beams = indices.size() + 1; // g0+
      A.resize(af::mat_grid(n_beams, n_beams));
      A(0, 0) = 0;
      for (size_t i = 1; i < n_beams; i++) {
        miller::index<> h_i = indices[i - 1];
        int ii = mi_lookup.find_hkl(h_i);
        SMTBX_ASSERT(ii >= 0);
        complex_t Fc_i = Fcs_k[ii];
        // h_i - (0,0,0)
        A(i, 0) = Fc_i;
        // (0,0,0) - h_i
        A(0, i) = std::conj(Fc_i);
        A(i, i) = 0;
        for (size_t j = i + 1; j < n_beams; j++) {
          miller::index<> h_j = indices[j - 1];
          int i_m_j = mi_lookup.find_hkl(h_i - h_j);
          SMTBX_ASSERT(i_m_j >= 0);
          complex_t Fc = Fcs_k[i_m_j];
          A(i, j) = Fc;
          A(j, i) = std::conj(Fc);
        }
      }
    }

    static bool sort_beams(const std::pair<size_t, FloatType>& a,
      const std::pair<size_t, FloatType>& b)
    {
      return (a.second < b.second);
    }

    /* Use_Sg controls how the beams are selected - using plain Sg or |Fc|/(Sg+wght)
    * In the case use_Sg is true - wght parameter is impose Sg limit on beams
    * added into the matrix
    */
    static af::shared<miller::index<> > build_Ug_matrix_N(cmat_t& A,
      const af::shared<complex_t>& Fcs_k,
      const lookup_t& mi_lookup,
      const af::shared<miller::index<> >& index_selection,
      FloatType Kl,
      const miller::index<> &h,
      const mat3_t &RMf,
      size_t num, bool use_Sg, FloatType wght)
    {
      typedef std::pair<size_t, FloatType> se_t;
      std::vector<se_t> beams;
      for (size_t i = 0; i < index_selection.size(); i++) {
        miller::index<> h_ = index_selection[i];
        if (h_ == h) {
          continue;
        }
        long idx = mi_lookup.find_hkl(h_);
        cart_t g = RMf * cart_t(h_[0], h_[1], h_[2]);
        FloatType Sg = std::abs(calc_Sg(g, Kl));
        if (use_Sg) {
          beams.push_back(std::make_pair(i, Sg));
        }
        else {
          beams.push_back(
            std::make_pair(i, std::abs(Fcs_k[idx]) / (Sg + wght)));
        }
      }
      std::sort(beams.begin(), beams.end(), sort_beams);
      af::shared<miller::index<> > indices;
      indices.push_back(h);
      if (use_Sg) {
        for (size_t i = 0; i < std::min(num, beams.size()); i++) {
          if (beams[i].second < 0.0025) {
            indices.push_back(index_selection[beams[i].first]);
            if (indices.size() == num) {
              break;
            }
          }
        }
      }
      else {
        for (size_t i = 0; i < std::min(num, beams.size()) - 1; i++) {
          indices.push_back(index_selection[beams[beams.size() - i - 1].first]);
        }
      }
      build_Ug_matrix(A, Fcs_k, mi_lookup, indices);
      return indices;
    }

    static void build_D_matrices(
      const lookup_t& mi_lookup,
      const af::shared<miller::index<> >& indices,
      const cmat_t& DM_kin,
      af::shared<cmat_t>& Ds_kin)
    {
      const size_t n_beams = indices.size() + 1; // g0+
      size_t n_param = DM_kin.accessor().n_columns();
      Ds_kin.reserve(n_param);
      for (size_t pi = 0; pi < n_param; pi++) {
        cmat_t D(af::mat_grid(n_beams, n_beams));
        for (size_t i = 1; i < n_beams; i++) {
          miller::index<> h_i = indices[i - 1];
          int ii = mi_lookup.find_hkl(h_i);
          SMTBX_ASSERT(ii >= 0);
          // h_i - (0,0,0)
          D(i, 0) = DM_kin(ii, pi);
          // (0,0,0) - h_i
          ii = mi_lookup.find_hkl(-h_i);
          SMTBX_ASSERT(ii >= 0);
          D(0, i) = DM_kin(ii, pi);
          for (size_t j = i + 1; j < n_beams; j++) {
            miller::index<> h_j = indices[j - 1];
            int i_m_j = mi_lookup.find_hkl(h_i - h_j);
            SMTBX_ASSERT(i_m_j >= 0);
            D(i, j) = DM_kin(i_m_j, pi);
            int j_m_i = mi_lookup.find_hkl(h_j - h_i);
            SMTBX_ASSERT(j_m_i >= 0);
            D(j, i) = DM_kin(j_m_i, pi);
          }
        }
        Ds_kin.push_back(D);
      }
    }

    // considers the original matrix Hermitian
    static void modify_Ug_matrix(
      cmat_t& A,
      af::shared<miller::index<> > const& s_indices, // matrix beams
      af::shared<complex_t> const& Fcs_k,
      lookup_t const& mi_lookup,
      af::shared<miller::index<> > const& w_indices, // weak beams for pertubation
      af::shared<FloatType> const& Exitations) // 2KSg for w_indices
    {
      const size_t n_beams = A.accessor().n_columns();
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
        A(i, i) -= d_mod;
        A(i, 0) -= u_mod;
      }
    }

    /* J. Appl. Cryst. (2022). 55 */
    static void build_eigen_matrix_recipro(
      cmat_t& A,
      af::shared<miller::index<> > const& indices,
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N,
      af::shared<FloatType>& Pgs)
    {
      // Projection of K onto normal of frame normal, K*frame.normal
      const FloatType Kn = N * K;
      const size_t n_beams = indices.size() + 1;
      A.resize(af::mat_grid(n_beams, n_beams));
      Pgs.resize(n_beams);
      for (size_t i = 0; i < n_beams; i++) {
        miller::index<> h_i = i == 0 ? miller::index<>(0, 0, 0) : indices[i - 1];
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
      cmat_t& A, // Ug matrix modified
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

    static af::shared<complex_t> calc_amps_recipro(
      cmat_t& A,
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
      cmat_t Bi = B.deep_copy();
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
      const complex_t exp_k(0, 2 * scitbx::constants::pi * thickness);
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
      cmat_t& A,
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
      cmat_t Bi = B.deep_copy();
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

    static complex_t calc_amp_2beam(
      const miller::index<>& h, complex_t Ug_,
      FloatType thickness,
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N,
      size_t idx=1)
    {
      FloatType Kn = N * K, Kl = K.length();
      cart_t K_g = K + RMf * cart_t(h[0], h[1], h[2]);
      FloatType K_gn = K_g * N;
      FloatType two_K_cos = 2 * Kl * K_gn / K_g.length();
      FloatType Sg = (Kl * Kl - K_g.length_sq()) / (2 * Kl);

      complex_t exp_k(0, 2 * scitbx::constants::pi * thickness);
      complex_t Ug = Ug_ / two_K_cos;
      complex_t A[4] = { 0,std::conj(Ug), Ug, Sg };
      FloatType v[2];

      //using namespace fast_linalg;
      //lapack_int info = heev(LAPACK_ROW_MAJOR, 'V', LAPACK_UPPER, 2,
      //  &A[0][0], 2, &ev[0]);
      //SMTBX_ASSERT(!info)(info);
      // A[1] and A[3] are real numbers
      math_utils<FloatType>::two_beam_eigen(&A[0], &v[0]);
      // incident beam
      if (idx == 0) {
        return A[0] * std::exp(v[0] * exp_k) * std::conj(A[0]) +
          A[1].real() * std::exp(v[1] * exp_k) * A[1].real();
      }
      return A[2] * std::exp(v[0] * exp_k) * std::conj(A[0]) +
        A[3].real() * std::exp(v[1] * exp_k) * A[1].real();
    }

    static FloatType calc_Sg(cart_t const& g,
      cart_t const& K)
    {
      FloatType Kl = K.length();
      cart_t Kg = K + g;
      return (Kl * Kl - Kg.length_sq()) / (2 * Kl);
    }

    // considers K being at (0, 0, -Kl)
    static FloatType calc_Sg(cart_t const& g,
      FloatType Kl)
    {
      FloatType Kg_sql = cart_t(g[0], g[1], g[2] - Kl).length_sq();
      return (Kl * Kl - Kg_sql) / (2 * Kl);
    }

    static bool is_excited_g(cart_t const& g,
      cart_t const& K, FloatType MaxSg, FloatType MaxG, FloatType precession_angle)
    {
      FloatType gl = g.length(),
        Sg = std::abs(calc_Sg(g, K));
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
    * The indices will fulfil the MaxSg and MaxG parameters and sorted by Sg
    */
    static af::shared<ExcitedBeam> generate_index_set(
      mat3_t const& RMf, // matrix to orthogonalise and rotate into the frame basis
      cart_t const& K,
      FloatType min_d,
      FloatType MaxG, FloatType MaxSg,
      uctbx::unit_cell const& unit_cell)
    {
      using namespace cctbx::miller;

      index_generator h_generator(unit_cell, sgtbx::space_group_type("P1"),
        true,
        min_d, false);

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
        //FloatType w = s * s * Kg_sq;
        all.push_back(ExcitedBeam(h, g, Sg, Sg));
      }
      std::sort(all.begin(), all.end(), &ExcitedBeam::compare);
      return all;
    }

    static af::shared<ExcitedBeam> generate_index_set(
      af::shared<mat3_t> const& RMfs, // matrix to orthogonalise and rotate into the frame basis
      cart_t const& K,
      FloatType min_d,
      FloatType MaxG, FloatType MaxSg,
      uctbx::unit_cell const& unit_cell)
    {
      using namespace cctbx::miller;

      index_generator h_generator(unit_cell, sgtbx::space_group_type("P1"),
        true,
        min_d, false);

      miller::index<> h;
      af::shared<ExcitedBeam> all;

      const FloatType max_f_sq = MaxG * MaxG,
        Kl = K.length(),
        Kl_sq = Kl * Kl;
      while (!(h = h_generator.next()).is_zero()) {
        FloatType min_sg = MaxSg;
        cart_t g;
        for (size_t mi = 0; mi < RMfs.size(); mi++) {
          g = RMfs[mi] * cart_t(h[0], h[1], h[2]);
          if (mi == 0) {
            FloatType g_sq = g.length_sq();
            if (g_sq > max_f_sq) {
              break;
            }
          }
          FloatType Kg_sq = (K + g).length_sq();
          FloatType s = Kl_sq - Kg_sq;
          FloatType Sg = std::abs(s / (2 * Kl));
          if (MaxSg > 0 && Sg > MaxSg) {
            continue;
          }
          if (Sg < min_sg) {
            min_sg = Sg;
          }
        }
        if (min_sg < MaxSg) {
          all.push_back(ExcitedBeam(h, g, min_sg, min_sg));
        }
      }
      std::sort(all.begin(), all.end(), &ExcitedBeam::compare);
      return all;
    }
  }; //struct smtbx::ED::utils
}} // namespace smtbx::ED