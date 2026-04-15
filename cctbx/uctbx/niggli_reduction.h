// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Derived from gemmi::GruberVector in include/gemmi/cellred.hpp
// Copyright 2017-2024 Global Phasing Ltd.  MPL-2.0.

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
   * The main loop mirrors the Python krivy_gruber_1976.reduction.step() function
   * exactly: A1 (N1 sort), A2 (N2 sort, restart if fired), A3/A4 (N3 sign,
   * always), A5-A8 (reduction, restart if fired).  n_iterations() counts
   * individual actions (cb_update calls) to match the Python implementation.
   *
   * The accumulated change-of-basis matrix r_inv is stored as a 3x3 integer
   * matrix (denominator = 1) in the same convention as fast_minimum_reduction
   * and reduction_base._r_inv.
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

        // Main reduction loop, matching Python's: while (self.step()): pass
        while (step()) {}
      }

      scitbx::mat3<int> const& r_inv() const { return r_inv_; }

      int n_iterations() const { return n_iterations_; }

      double epsilon() const { return epsilon_; }

      //! Direct access to the reduced G6 Gruber parameters.
      double A()    const { return A_; }
      double B()    const { return B_; }
      double C()    const { return C_; }
      double xi()   const { return xi_; }
      double eta()  const { return eta_; }
      double zeta() const { return zeta_; }

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

      // --- Iteration counter (mirrors Python cb_update) ---
      // Called once per action, checks limit before incrementing.

      void cb_update()
      {
        if (n_iterations_ == iteration_limit_) {
          throw error_niggli_iteration_limit_exceeded();
        }
        n_iterations_++;
      }

      // --- Individual action methods (A1-A8 = N1, N2, N3, A5-A8) ---
      // Each mirrors the corresponding Python action method.

      // A1: N1 sort (A <= B).  Returns true if the swap was applied.
      bool do_a1()
      {
        if (eps_gt(A_, B_) ||
            (eps_eq(A_, B_) && eps_gt(std::abs(xi_), std::abs(eta_)))) {
          std::swap(A_, B_);
          std::swap(xi_, eta_);
          swap_columns_and_negate(0, 1);
          cb_update();
          return true;
        }
        return false;
      }

      // A2: N2 sort (B <= C).  Returns true if the swap was applied.
      bool do_a2()
      {
        if (eps_gt(B_, C_) ||
            (eps_eq(B_, C_) && eps_gt(std::abs(eta_), std::abs(zeta_)))) {
          std::swap(B_, C_);
          std::swap(eta_, zeta_);
          swap_columns_and_negate(1, 2);
          cb_update();
          return true;
        }
        return false;
      }

      // Helper: classify xi_/eta_/zeta_ signs.
      // n_positive: number strictly greater than +epsilon
      // n_nonneg:   number greater than -epsilon (i.e., >= 0 within tolerance)
      void sign_counts(int& n_positive, int& n_nonneg) const
      {
        n_positive = (xi_   >  epsilon_) + (eta_  >  epsilon_) + (zeta_ >  epsilon_);
        n_nonneg   = (xi_   > -epsilon_) + (eta_  > -epsilon_) + (zeta_ > -epsilon_);
      }

      // A3/A4: N3 sign normalization.  Always applied (always calls cb_update).
      // Mirrors Python n3_true_action (def_gt_0=true) and n3_false_action (false).
      void do_a3_a4()
      {
        int pos, nonneg;
        sign_counts(pos, nonneg);

        // Target sign: +1 if all-positive case (def_gt_0), -1 otherwise.
        // def_gt_0 <==> n_positive==3 OR (n_zero==0 AND n_positive==1)
        //           <==> pos==nonneg AND pos%2==1
        double sgn = (pos == nonneg && pos % 2 == 1) ? 1.0 : -1.0;

        // Negate columns whose current sign opposes the target.
        if (sgn * xi_   < -epsilon_) negate_column(0);
        if (sgn * eta_  < -epsilon_) negate_column(1);
        if (sgn * zeta_ < -epsilon_) negate_column(2);

        // Edge case: zero-valued component in the all-negative branch requires
        // one extra negation to maintain the correct sign product
        // (mirrors Python n3_false_action's f[z]=-1 adjustment).
        if (pos != nonneg && pos % 2 == 1) {
          negate_column(
            std::fabs(zeta_) <= epsilon_ ? 2 :
            std::fabs(eta_)  <= epsilon_ ? 1 : 0);
        }

        xi_   = std::copysign(xi_,   sgn);
        eta_  = std::copysign(eta_,  sgn);
        zeta_ = std::copysign(zeta_, sgn);

        cb_update();  // ONE action regardless of how many columns were negated
      }

      // A5: reduce |xi| w.r.t. B.  Returns true if the step was applied.
      bool do_a5()
      {
        if (eps_gt(std::abs(xi_), B_) ||
            (xi_ >= B_ - epsilon_ && eps_lt(2*eta_, zeta_)) ||
            (xi_ <= -(B_ - epsilon_) && eps_lt(zeta_, 0.0))) {
          double s = xi_ >= 0 ? 1.0 : -1.0;
          C_   += B_ - xi_ * s;
          eta_ -= zeta_ * s;
          xi_  -= 2 * B_ * s;
          add_column(1, 2, -(int)s);
          cb_update();
          return true;
        }
        return false;
      }

      // A6: reduce |eta| w.r.t. A.  Returns true if the step was applied.
      bool do_a6()
      {
        if (eps_gt(std::abs(eta_), A_) ||
            (eta_ >= A_ - epsilon_ && eps_lt(2*xi_, zeta_)) ||
            (eta_ <= -(A_ - epsilon_) && eps_lt(zeta_, 0.0))) {
          double s = eta_ >= 0 ? 1.0 : -1.0;
          C_   += A_ - eta_ * s;
          xi_  -= zeta_ * s;
          eta_ -= 2 * A_ * s;
          add_column(0, 2, -(int)s);
          cb_update();
          return true;
        }
        return false;
      }

      // A7: reduce |zeta| w.r.t. A.  Returns true if the step was applied.
      bool do_a7()
      {
        if (eps_gt(std::abs(zeta_), A_) ||
            (zeta_ >= A_ - epsilon_ && eps_lt(2*xi_, eta_)) ||
            (zeta_ <= -(A_ - epsilon_) && eps_lt(eta_, 0.0))) {
          double s = zeta_ >= 0 ? 1.0 : -1.0;
          B_    += A_ - zeta_ * s;
          xi_   -= eta_ * s;
          zeta_ -= 2 * A_ * s;
          add_column(0, 1, -(int)s);
          cb_update();
          return true;
        }
        return false;
      }

      // A8: mixed condition.  Returns true if the step was applied.
      bool do_a8()
      {
        if (eps_lt(xi_ + eta_ + zeta_ + A_ + B_, 0.0) ||
            (eps_eq(xi_ + eta_ + zeta_ + A_ + B_, 0.0) &&
             eps_gt(2*(A_ + eta_) + zeta_, 0.0))) {
          C_   += A_ + B_ + xi_ + eta_ + zeta_;
          xi_  += 2*B_ + zeta_;
          eta_ += 2*A_ + zeta_;
          add_column(0, 2, 1);
          add_column(1, 2, 1);
          cb_update();  // ONE action for the two add_column calls (mirrors Python)
          return true;
        }
        return false;
      }

      // --- Main step function ---
      // Mirrors Python krivy_gruber_1976.reduction.step() exactly:
      //   A1 (optional N1 sort)
      //   A2 (N2 sort, restart if fired)
      //   A3/A4 (N3 sign, always)
      //   A5-A8 (reduction, restart if fired)
      //   return false when nothing triggered a restart

      bool step()
      {
        do_a1();                    // A1: N1 sort (optional, no restart)
        if (do_a2()) return true;   // A2: N2 sort (restart if fired)
        do_a3_a4();                 // A3/A4: sign normalization (always)
        if (do_a5()) return true;   // A5
        if (do_a6()) return true;   // A6
        if (do_a7()) return true;   // A7
        if (do_a8()) return true;   // A8
        return false;               // done
      }

      double A_, B_, C_, xi_, eta_, zeta_;
      double epsilon_;
      int n_iterations_;
      int iteration_limit_;
      scitbx::mat3<int> r_inv_;
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_NIGGLI_REDUCTION_H
