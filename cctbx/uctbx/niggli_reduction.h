#ifndef CCTBX_UCTBX_NIGGLI_REDUCTION_H
#define CCTBX_UCTBX_NIGGLI_REDUCTION_H

#include <cctbx/uctbx.h>
#include <cctbx/error.h>
#include <scitbx/mat3.h>
#include <cmath>
#include <algorithm>

namespace cctbx { namespace uctbx {

  class error_niggli_iteration_limit_exceeded : public error
  {
    public:
      error_niggli_iteration_limit_exceeded()
      :
        error("Niggli reduction iteration limit exceeded.")
      {}
  };

  //! C++ implementation of the Krivy-Gruber (1976) Niggli reduction.
  /*!
   * Implements the algorithm of Krivy & Gruber (1976) Acta Cryst. A32, 297-298.
   * Operates on the G6 Gruber parameterization:
   *   A = a.a,  B = b.b,  C = c.c
   *   xi = 2(b.c),  eta = 2(a.c),  zeta = 2(a.b)
   *
   * The accumulated change-of-basis matrix r_inv is stored as a 3x3 integer
   * matrix (denominator = 1, entries are 0 or +/-1) in the same convention as
   * fast_minimum_reduction::r_inv() and reduction_base._r_inv.
   *
   * Based on the algorithm in gemmi::GruberVector (cellred.hpp), ported to
   * cctbx to remove the gemmi external dependency.
   */
  class niggli_reduction
  {
    public:
      niggli_reduction(
        uctbx::unit_cell const& unit_cell,
        double relative_epsilon=1.e-5,
        int iteration_limit=1000)
      :
        n_iterations_(0),
        iteration_limit_(iteration_limit),
        r_inv_(1,0,0, 0,1,0, 0,0,1)
      {
        // Initialize G6 from the metrical matrix.
        // metrical_matrix() returns sym_mat3 with indices:
        //   [0]=a*a, [1]=b*b, [2]=c*c, [3]=a*b, [4]=a*c, [5]=b*c
        scitbx::sym_mat3<double> const& m = unit_cell.metrical_matrix();
        A_    = m[0];
        B_    = m[1];
        C_    = m[2];
        xi_   = 2.0 * m[5];  // 2*(b.c)
        eta_  = 2.0 * m[4];  // 2*(a.c)
        zeta_ = 2.0 * m[3];  // 2*(a.b)
        epsilon_ = std::pow(unit_cell.volume(), 1.0/3.0) * relative_epsilon;

        // Main reduction loop: normalize then apply one Niggli step.
        for (;;) {
          normalize();
          if (++n_iterations_ == iteration_limit_) {
            throw error_niggli_iteration_limit_exceeded();
          }
          if (niggli_step()) {
            break;
          }
        }
      }

      scitbx::mat3<int> const& r_inv() const { return r_inv_; }

      int n_iterations() const { return n_iterations_; }

      uctbx::unit_cell as_unit_cell() const
      {
        // as_sym_mat3 convention: (A, B, C, zeta/2, eta/2, xi/2)
        // which maps to metrical matrix (a*a, b*b, c*c, a*b, a*c, b*c)
        scitbx::sym_mat3<double> mm(A_, B_, C_, zeta_/2, eta_/2, xi_/2);
        return uctbx::unit_cell(mm);
      }

    protected:

      // --- Column operations on r_inv_ ---
      // Columns correspond to basis vectors; all reduction steps are
      // column operations on the change-of-basis matrix.

      void swap_columns_and_negate(int ci, int cj)
      {
        for (int row = 0; row < 3; row++) {
          std::swap(r_inv_[row*3+ci], r_inv_[row*3+cj]);
        }
        for (int k = 0; k < 9; k++) {
          r_inv_[k] = -r_inv_[k];
        }
      }

      void negate_column(int ci)
      {
        for (int row = 0; row < 3; row++) {
          r_inv_[row*3+ci] = -r_inv_[row*3+ci];
        }
      }

      void add_column(int src, int dst, int sign)
      {
        for (int row = 0; row < 3; row++) {
          r_inv_[row*3+dst] += sign * r_inv_[row*3+src];
        }
      }

      // --- Epsilon-based comparisons ---

      bool eps_lt(double x, double y) const { return x < y - epsilon_; }
      bool eps_gt(double x, double y) const { return eps_lt(y, x); }
      bool eps_eq(double x, double y) const
      {
        return !eps_lt(x, y) && !eps_lt(y, x);
      }

      // --- N1 sorting step (used inside normalize) ---

      void do_step_N1()
      {
        if (eps_gt(A_, B_) ||
            (eps_eq(A_, B_) && eps_gt(std::abs(xi_), std::abs(eta_)))) {
          std::swap(A_, B_);
          std::swap(xi_, eta_);
          swap_columns_and_negate(0, 1);
        }
      }

      // --- normalize(): Algorithm N from Gruber (1973) ---
      // Sorts A <= B <= C with tie-breaking on |xi|, |eta|, |zeta|,
      // then enforces consistent signs on xi, eta, zeta.

      void normalize()
      {
        // N1: sort A <= B
        do_step_N1();

        // N2: sort B <= C (then re-apply N1 since B may have decreased)
        if (eps_gt(B_, C_) ||
            (eps_eq(B_, C_) && eps_gt(std::abs(eta_), std::abs(zeta_)))) {
          std::swap(B_, C_);
          std::swap(eta_, zeta_);
          swap_columns_and_negate(1, 2);
          do_step_N1();
        }

        // N3: enforce consistent sign convention for xi, eta, zeta.
        // Target: all three values have the same sign.
        // sgn=+1 if exactly an odd number are strictly positive (all zeros ignored),
        //        AND there are no zero-valued components.
        // sgn=-1 otherwise (make all non-positive).
        int pos   = (xi_ > epsilon_)    + (eta_ > epsilon_)    + (zeta_ > epsilon_);
        int nonneg = (xi_ >= -epsilon_) + (eta_ >= -epsilon_) + (zeta_ >= -epsilon_);
        double sgn = (pos == nonneg && pos % 2 == 1) ? 1.0 : -1.0;

        // Negate columns whose value has the wrong sign.
        if (sgn * xi_   < -epsilon_) negate_column(0);
        if (sgn * eta_  < -epsilon_) negate_column(1);
        if (sgn * zeta_ < -epsilon_) negate_column(2);

        // Special case: if there is a zero component and we are in the odd-positive
        // case, negate one more column to keep the product of all signs consistent.
        if (pos != nonneg && pos % 2 == 1) {
          negate_column(
            std::fabs(zeta_) <= epsilon_ ? 2 :
            std::fabs(eta_)  <= epsilon_ ? 1 : 0);
        }

        // Apply the target sign to the scalar values.
        xi_   = std::copysign(xi_,   sgn);
        eta_  = std::copysign(eta_,  sgn);
        zeta_ = std::copysign(zeta_, sgn);
      }

      // --- niggli_step(): Steps A5-A8 from Krivy & Gruber (1976) ---
      // Must be called after normalize(). Returns true if the cell is already
      // Niggli-reduced (no step was needed).

      bool niggli_step()
      {
        // Step 5: reduce |xi| w.r.t. B
        if (eps_gt(std::abs(xi_), B_) ||
            (xi_ >= B_ - epsilon_ && eps_lt(2*eta_, zeta_)) ||
            (xi_ <= -(B_ - epsilon_) && eps_lt(zeta_, 0.0))) {
          double s = xi_ >= 0 ? 1.0 : -1.0;
          C_   += B_ - xi_ * s;
          eta_ -= zeta_ * s;
          xi_  -= 2 * B_ * s;
          add_column(1, 2, -(int)s);
          return false;
        }
        // Step 6: reduce |eta| w.r.t. A
        if (eps_gt(std::abs(eta_), A_) ||
            (eta_ >= A_ - epsilon_ && eps_lt(2*xi_, zeta_)) ||
            (eta_ <= -(A_ - epsilon_) && eps_lt(zeta_, 0.0))) {
          double s = eta_ >= 0 ? 1.0 : -1.0;
          C_   += A_ - eta_ * s;
          xi_  -= zeta_ * s;
          eta_ -= 2 * A_ * s;
          add_column(0, 2, -(int)s);
          return false;
        }
        // Step 7: reduce |zeta| w.r.t. A
        if (eps_gt(std::abs(zeta_), A_) ||
            (zeta_ >= A_ - epsilon_ && eps_lt(2*xi_, eta_)) ||
            (zeta_ <= -(A_ - epsilon_) && eps_lt(eta_, 0.0))) {
          double s = zeta_ >= 0 ? 1.0 : -1.0;
          B_    += A_ - zeta_ * s;
          xi_   -= eta_ * s;
          zeta_ -= 2 * A_ * s;
          add_column(0, 1, -(int)s);
          return false;
        }
        // Step 8: mixed condition
        if (eps_lt(xi_ + eta_ + zeta_ + A_ + B_, 0.0) ||
            (eps_eq(xi_ + eta_ + zeta_ + A_ + B_, 0.0) &&
             eps_gt(2*(A_ + eta_) + zeta_, 0.0))) {
          C_   += A_ + B_ + xi_ + eta_ + zeta_;
          xi_  += 2*B_ + zeta_;
          eta_ += 2*A_ + zeta_;
          add_column(0, 2, 1);
          add_column(1, 2, 1);
          return false;
        }
        return true;  // already Niggli-reduced
      }

      double A_, B_, C_, xi_, eta_, zeta_;
      double epsilon_;
      int n_iterations_;
      int iteration_limit_;
      scitbx::mat3<int> r_inv_;
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_NIGGLI_REDUCTION_H
