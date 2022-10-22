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

    static void build_eigen_matrix(
      af::shared<complex_t> const& Fcs_k,
      lookup_t const& mi_lookup,
      af::shared<miller::index<> > indices, // implicit {0,0,0} at 0
      cart_t const& K,
      mat3_t const& RMf,
      cart_t const& N,
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<FloatType>& M,
      FloatType Fc2Ug)
    {
      // Projection of K onto normal of frame normal, K*frame.normal
      const FloatType Kn = N * K,
        Kl = K.length();
      using namespace fast_linalg;
      const size_t n_beams = indices.size() + 1; // g0+
      A.resize(af::mat_grid(n_beams, n_beams));
      M.resize(n_beams);
      M[0] = 1; //for g0
      A(0, 0) = 0;
      for (size_t i = 1; i < n_beams; i++) {
        miller::index<> h_i = indices[i - 1];
        int ii = mi_lookup.find_hkl(h_i);
        complex_t Fc_i = ii != -1 ? Fcs_k[ii] : 0;

        cart_t g_i = RMf * cart_t(h_i[0], h_i[1], h_i[2]);
        FloatType s = (Kl * Kl - (K + g_i).length_sq());
        FloatType i_den = std::sqrt(1. / (1 + g_i * N / Kn));
        A(i, i) = s * i_den * i_den;

        A(i, 0) = Fc2Ug * Fc_i * i_den;
        A(0, i) = std::conj(A(i, 0));

        M[i] = i_den;
        for (size_t j = i + 1; j < n_beams; j++) {
          miller::index<> h_j = indices[j - 1];
          cart_t g_j = RMf * cart_t(h_j[0], h_j[1], h_j[2]);
          miller::index<> h_i_m_j = h_i - h_j;
          int i_m_j = mi_lookup.find_hkl(h_i_m_j);
          complex_t Fc_i_m_j = i_m_j < 0 ? 0 : Fc_i_m_j = Fcs_k[i_m_j];
          FloatType j_den = std::sqrt(1. / (1 + g_j * N / Kn));
          A(i, j) = Fc2Ug * Fc_i_m_j * i_den * j_den;
          A(j, i) = std::conj(A(i, j));
        }
      }
    }

    /*
        // eigenvalues, this is for generic matrix
        af::shared<complex_t> ev(n_beams);
        // right eigenvectors
        af::versa<complex_t, af::c_grid<2> > eV(af::c_grid<2>(n_beams, n_beams));
        lapack_int info = geev(LAPACK_ROW_MAJOR, 'N', 'V', n_beams,
          A.begin(), n_beams, ev.begin(), 0, n_beams, eV.begin(), n_beams);
        SMTBX_ASSERT(!info)(info);
        af::versa<complex_t, af::mat_grid>
          B(af::mat_grid(n_beams, n_beams), *eV.begin());
        // invert eV now
        {
          af::shared<lapack_int> pivots(n_beams);
          info = getrf(LAPACK_ROW_MAJOR, n_beams, n_beams, eV.begin(),
            n_beams, pivots.begin());
          SMTBX_ASSERT(!info)(info);
          info = getri(LAPACK_ROW_MAJOR, n_beams, eV.begin(),
            n_beams, pivots.begin());
          SMTBX_ASSERT(!info)(info);
        }
    */

    static af::shared<FloatType> calc_I(
      af::versa<complex_t, af::mat_grid>& A,
      af::shared<FloatType> const& M,
      FloatType thickness,
      FloatType Kn,
      size_t num)
    {
      using namespace fast_linalg;
      const size_t n_beams = A.accessor().n_columns();
      af::shared<FloatType> ev(n_beams);
      lapack_int info = heev(LAPACK_ROW_MAJOR, 'V', LAPACK_UPPER, n_beams,
        A.begin(), n_beams, ev.begin());
      SMTBX_ASSERT(!info)(info);
      // heev replaces A with column-wise eigenvectors
      af::versa<complex_t, af::mat_grid> B(af::mat_grid(n_beams, n_beams));
      af::versa<complex_t, af::mat_grid> Bi(af::mat_grid(n_beams, n_beams));

      const complex_t exp_k(0, scitbx::constants::pi * thickness / Kn);
      // init Bi = B^H and update B = B*diag(exp(2*pi*thickness*i*ev/(2*Kn))
      for (size_t i = 0; i < n_beams; i++) {
        complex_t m = std::exp(ev[i] * exp_k);
        for (size_t j = 0; j < n_beams; j++) {
          B(j, i) = A(j, i) * m;
          Bi(i, j) = std::conj(A(j, i));
        }
      }
      // S = B*Bi
      af::versa<complex_t, af::mat_grid> S = af::matrix_multiply(B.const_ref(), Bi.const_ref());
      //update S = M*S*M^-1, m - diagonal, M[0] = 1, so start with 1
      for (size_t i = 1; i < n_beams; i++) {
        for (size_t j = 1; j < n_beams; j++) {
          S(i, j) *= M[i]; // M is on the left - apply to rows
          S(j, i) /= M[i]; // M^-1 is on the right - > apply to cols
        }
      }
      // extract first col of S
      //af::shared<complex_t> u(n_beams);
      //u(0) = 1;
      //af::shared<complex_t> up = af::matrix_multiply(S.const_ref(), u.const_ref());
      af::shared<FloatType> rv;
      for (size_t i = 1; i < n_beams; i++) {
        rv.push_back(std::norm(S(i, 0)));
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
      return Sg < MaxSg&& gl < MaxG&& Sg / (gl * precession_angle) < 1.0;
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
    * The indices are sorted by excitation FoM desc
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

      index<> h;
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