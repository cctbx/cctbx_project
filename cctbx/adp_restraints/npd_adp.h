#ifndef CCTBX_ADP_RESTRAINTS_NPD_ADP_H
#define CCTBX_ADP_RESTRAINTS_NPD_ADP_H

#include <cctbx/adp_restraints/adp_restraints.h>

namespace cctbx { namespace adp_restraints {

  /** A restraint that keeps an anisotropic ADP positive-definite, the
      SHELXL XNPD instruction.

      The physical quantity is the smallest eigenvalue of U, which is negative
      exactly when the ellipsoid is non-positive-definite (NPD). Rather than
      decompose U each cycle, this restrains a smooth, scale-invariant proxy
      for it:

          s(U) = det(U) / U_eq^3,   U_eq = trace(U)/3

      s = 1 for a sphere, s in (0, 1] for a positive-definite ellipsoid, s -> 0
      as one axis collapses, s < 0 once the atom goes NPD. Normalising by U_eq^3
      makes s dimensionless, so a single target means the same thing for a tight
      light atom and a large heavy one, which a bare det(U) (units Angstrom^6,
      scaling as size^3) would not.

      One-sided: the penalty w*(s_target - s)^2 is applied only to atoms
      currently below s_target, re-evaluated every cycle, so healthy atoms are
      untouched and only those heading NPD are pulled back. det(U) is smooth, so
      the only non-smoothness is the switch at the threshold, handled by adding
      the equation with zero weight when inactive.

      The target and the weight are carried on the proxy, so a caller can tune
      them.
   */
  struct npd_adp_proxy : public adp_restraint_proxy<1> {
    npd_adp_proxy() {}

    npd_adp_proxy(af::tiny<unsigned, 1> const &i_seq,
      double weight, double s_target)
    : adp_restraint_proxy<1>(i_seq, weight),
      s_target(s_target)
    {}

    //! Target for det(U)/U_eq^3, the margin above the NPD boundary at s = 0.
    double s_target;
  };

  class npd_adp {
  public:
    npd_adp(scitbx::sym_mat3<double> const &u_cart,
      double weight, double s_target)
    : weight(weight), s_target(s_target), use_u_aniso(true)
    {
      init(u_cart);
    }

    npd_adp(
      adp_restraint_params<double> const &params,
      npd_adp_proxy const &proxy)
    : weight(proxy.weight), s_target(proxy.s_target)
    {
      std::size_t i_seq = proxy.i_seqs[0];
      CCTBX_ASSERT(i_seq < params.use_u_aniso.size());
      use_u_aniso = params.use_u_aniso[i_seq];
      if (use_u_aniso) {
        CCTBX_ASSERT(i_seq < params.u_cart.size());
        init(params.u_cart[i_seq]);
      }
      else {
        // an isotropic atom cannot be NPD; nothing to restrain
        active_ = false;
        delta_ = 0;
        s_ = 1;
      }
    }

    //! Whether this atom is currently below the target and being restrained.
    bool active() const { return active_; }
    //! The proxy quantity det(U)/U_eq^3.
    double s() const { return s_; }
    //! s_target - s where restrained, else 0.
    double delta() const { return delta_; }
    //! weight * delta^2, zero unless restrained.
    double residual() const { return active_ ? weight * delta_ * delta_ : 0; }
    //! |delta|, for the sorted report.
    double rms_deltas() const { return std::abs(delta_); }

    //! d(residual)/d(u_cart) as a symmetric tensor, for the gradient-vector
    //! path (show_sorted, non-least-squares optimisers). The least-squares
    //! path uses linearise() below.
    scitbx::sym_mat3<double> gradients() const {
      scitbx::sym_mat3<double> g(0,0,0,0,0,0);
      if (!active_) return g;
      // residual = w delta^2, d/du = 2 w delta * d(delta)/du = 2 w delta * cart_grad_
      for (int i=0; i<6; i++) g[i] = 2 * weight * delta_ * cart_grad_[i];
      return g;
    }

    void add_gradients(
      af::ref<scitbx::sym_mat3<double> > const &gradients_aniso_cart,
      af::tiny<unsigned, 1> const &i_seqs) const
    {
      if (active_) gradients_aniso_cart[i_seqs[0]] += gradients();
    }

    void linearise(
      uctbx::unit_cell const &unit_cell,
      cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
      cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const
        &parameter_map,
      af::tiny<unsigned, 1> const &i_seqs)
    {
      // One row per proxy either way, so the row total the manager pre-counted
      // is exact; an inactive atom contributes a zero-weight row, which adds
      // nothing to the normal equations.
      std::size_t row_i = linearised_eqns.next_row();
      if (!active_) {
        linearised_eqns.weights[row_i] = 0;
        linearised_eqns.deltas[row_i] = 0;
        return;
      }
      cctbx::xray::parameter_indices const &ids = parameter_map[i_seqs[0]];
      CCTBX_ASSERT(ids.u_aniso != -1);
      // delta is a function of u_cart; the refined parameters are u_star, so
      // d(delta)/d(u_star) = M^T d(delta)/d(u_cart), and then the off-diagonal
      // columns carry a factor of two. That factor is the cctbx-wide convention
      // for how a u_star gradient enters the design matrix -- linearise_1 does
      // exactly this, and the finite-difference check in tst_npd_adp.py doubles
      // its numerical off-diagonal gradient to match, so it is not optional.
      scitbx::sym_mat3<double> grad_u_star;
      scitbx::matrix::matrix_transposed_vector(
        6, 6, unit_cell.u_star_to_u_cart_linear_map().begin(),
        cart_grad_.begin(), grad_u_star.begin());
      for (int j=0; j<6; j++) {
        linearised_eqns.design_matrix(row_i, ids.u_aniso+j) =
          (j > 2 ? 2*grad_u_star[j] : grad_u_star[j]);
      }
      linearised_eqns.weights[row_i] = weight;
      linearised_eqns.deltas[row_i] = delta_;
    }

    double weight;
    double s_target;
    bool use_u_aniso;

  protected:
    void init(scitbx::sym_mat3<double> const &u) {
      double const a = u[0], b = u[1], c = u[2], d = u[3], e = u[4], f = u[5];
      double const u_eq = (a + b + c) / 3;
      // U_eq is positive for any physical atom, even a slightly NPD one whose
      // trace is still positive; guard the division for a grossly wrong atom
      // rather than divide by zero.
      if (u_eq <= 1e-12) {
        active_ = false;
        delta_ = 0;
        s_ = -1;   // definitely not positive definite
        return;
      }
      double const det = a*b*c + 2*d*e*f - a*f*f - b*e*e - c*d*d;
      double const ueq3 = u_eq*u_eq*u_eq;
      s_ = det / ueq3;
      active_ = use_u_aniso && (s_ < s_target);
      if (!active_) {
        delta_ = 0;
        return;
      }
      delta_ = s_target - s_;
      // d(det)/d(component), symmetric-matrix determinant differentiated with
      // respect to each of the six independent components (the off-diagonal
      // factor of two falls out of the det expression itself):
      double const det_a = b*c - f*f;
      double const det_b = a*c - e*e;
      double const det_c = a*b - d*d;
      double const det_d = 2*(e*f - c*d);
      double const det_e = 2*(d*f - b*e);
      double const det_f = 2*(d*e - a*f);
      // s = det * u_eq^-3; d(s)/dk = det_k/u_eq^3 - 3 s u_eq^-1 d(u_eq)/dk,
      // with d(u_eq)/dk = 1/3 for the diagonal components and 0 otherwise.
      double const inv_ueq3 = 1 / ueq3;
      double const s_over_ueq = s_ / u_eq;
      scitbx::sym_mat3<double> ds;
      ds[0] = det_a*inv_ueq3 - s_over_ueq;
      ds[1] = det_b*inv_ueq3 - s_over_ueq;
      ds[2] = det_c*inv_ueq3 - s_over_ueq;
      ds[3] = det_d*inv_ueq3;
      ds[4] = det_e*inv_ueq3;
      ds[5] = det_f*inv_ueq3;
      // delta = s_target - s, so d(delta)/d(u_cart) = -d(s)/d(u_cart)
      for (int i=0; i<6; i++) cart_grad_[i] = -ds[i];
    }

    bool active_;
    double delta_;
    double s_;
    scitbx::sym_mat3<double> cart_grad_;  // d(delta)/d(u_cart), convention A
  };

}} // cctbx::adp_restraints

#endif // GUARD
