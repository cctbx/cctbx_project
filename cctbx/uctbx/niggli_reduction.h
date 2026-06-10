#ifndef CCTBX_UCTBX_NIGGLI_REDUCTION_H
#define CCTBX_UCTBX_NIGGLI_REDUCTION_H

#include <cctbx/uctbx.h>
#include <cctbx/error.h>
#include <scitbx/array_family/tiny.h>
#include <cmath>
#include <string>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/sgtbx/rot_mx.h>

namespace cctbx { namespace uctbx {

  //! Exception raised when the Niggli reduction iteration limit is exceeded.
  class niggli_reduction_iteration_limit_exceeded : public error
  {
    public:
      niggli_reduction_iteration_limit_exceeded(std::size_t limit)
      :
        error(
          "Niggli reduction iteration limit exceeded (limit="
          + std::to_string(limit) + ").")
      {}
  };

  /*! \brief C++ implementation of the Krivy-Gruber 1976 Niggli cell
      reduction algorithm.

      Reference: I. Krivy and B. Gruber, Acta Cryst. A32, 297-298 (1976).

      Uses epsilon-aware comparisons throughout (unlike fast_minimum_reduction
      which uses exact comparisons).  The algorithm and action matrices are an
      exact port of cctbx/uctbx/krivy_gruber_1976.py.

      Template parameters:
        FloatType       – floating-point type for cell parameters (default double)
        IntFromFloatType – integer type for change-of-basis matrix (default int)
   */
  template <typename FloatType=double,
            typename IntFromFloatType=int>
  class niggli_reduction
  {
    public:
      /*! \brief Execute the Niggli reduction algorithm.

          \param unit_cell      Input cell to reduce.
          \param relative_epsilon  Tolerance relative to volume^(1/3).
                                   Caller must resolve Python None to 1.e-5.
          \param iteration_limit   Maximum number of change-of-basis steps.
                                   Caller must resolve Python None to 1000.
       */
      niggli_reduction(
        uctbx::unit_cell const& unit_cell,
        FloatType relative_epsilon,
        std::size_t iteration_limit)
      :
        iteration_limit_(iteration_limit),
        r_inv_(scitbx::mat3<IntFromFloatType>(1)),
        n_iterations_(0)
      {
        uc_sym_mat3 const& g = unit_cell.metrical_matrix();
        a_ = g[0];
        b_ = g[1];
        c_ = g[2];
        d_ = 2 * g[5];
        e_ = 2 * g[4];
        f_ = 2 * g[3];
        FloatType vol = unit_cell.volume();
        epsilon_ = std::pow(vol, FloatType(1) / FloatType(3))
                   * relative_epsilon;
        while (step());
      }

      //! Gruber parameterization (a, b, c, d, e, f).
      af::double6
      as_gruber_matrix() const
      {
        return af::double6(a_, b_, c_, d_, e_, f_);
      }

      //! Niggli parameterization (a, b, c, d/2, e/2, f/2).
      af::double6
      as_niggli_matrix() const
      {
        return af::double6(a_, b_, c_, d_/2, e_/2, f_/2);
      }

      //! Parameterization compatible with uctbx::unit_cell::metrical_matrix().
      uc_sym_mat3
      as_sym_mat3() const
      {
        return uc_sym_mat3(a_, b_, c_, f_/2, e_/2, d_/2);
      }

      //! Conversion to uctbx::unit_cell.
      uctbx::unit_cell
      as_unit_cell() const
      {
        return uctbx::unit_cell(as_sym_mat3());
      }

      //! Change-of-basis matrix compatible with uctbx::unit_cell::change_basis().
      scitbx::mat3<IntFromFloatType> const&
      r_inv() const { return r_inv_; }

      //! Number of change-of-basis steps executed.
      std::size_t
      n_iterations() const { return n_iterations_; }

      //! Value as passed to the constructor.
      std::size_t
      iteration_limit() const { return iteration_limit_; }

      //! Cell type (1 = all d,e,f positive; 2 = all non-positive; 0 = mixed).
      int
      type() const
      {
        af::tiny<int, 2> nz_np = def_test();
        if (nz_np[1] == 3) return 1;
        if (nz_np[1] == 0) return 2;
        return 0;
      }

      //! True if a <= b <= c and |d|<=b, |e|<=a, |f|<=a (epsilon-aware).
      bool
      meets_primary_conditions() const
      {
        if (eps_gt(a_, b_)) return false;
        if (eps_gt(b_, c_)) return false;
        if (eps_gt(std::abs(d_), b_)) return false;
        if (eps_gt(std::abs(e_), a_)) return false;
        if (eps_gt(std::abs(f_), a_)) return false;
        return true;
      }

      //! True if primary conditions hold, type != 0, and type-2 sum condition.
      bool
      meets_main_conditions() const
      {
        if (!meets_primary_conditions()) return false;
        int t = type();
        if (t == 0) return false;
        if (t == 2) {
          if (eps_lt(d_ + e_ + f_ + a_ + b_, FloatType(0))) return false;
        }
        return true;
      }

      //! True if meets_main_conditions and tie-breaking on equal lengths.
      bool
      is_buerger_cell() const
      {
        if (!meets_main_conditions()) return false;
        if (eps_eq(a_, b_)) {
          if (eps_gt(std::abs(d_), std::abs(e_))) return false;
        }
        if (eps_eq(b_, c_)) {
          if (eps_gt(std::abs(e_), std::abs(f_))) return false;
        }
        return true;
      }

      //! Change-of-basis operator that transforms the input cell to the Niggli cell.
      sgtbx::change_of_basis_op
      change_of_basis_op() const
      {
        return sgtbx::change_of_basis_op(
          sgtbx::rt_mx(sgtbx::rot_mx(r_inv_, 1))).inverse();
      }

      //! True if Buerger cell and 7 Niggli boundary uniqueness conditions.
      bool
      is_niggli_cell() const
      {
        if (!is_buerger_cell()) return false;
        if (eps_eq(d_, b_)) {
          if (eps_gt(f_, e_ + e_)) return false;
        }
        if (eps_eq(e_, a_)) {
          if (eps_gt(f_, d_ + d_)) return false;
        }
        if (eps_eq(f_, a_)) {
          if (eps_gt(e_, d_ + d_)) return false;
        }
        if (eps_eq(d_, -b_)) {
          if (!eps_eq(f_, FloatType(0))) return false;
        }
        if (eps_eq(e_, -a_)) {
          if (!eps_eq(f_, FloatType(0))) return false;
        }
        if (eps_eq(f_, -a_)) {
          if (!eps_eq(e_, FloatType(0))) return false;
        }
        if (eps_eq(d_ + e_ + f_ + a_ + b_, FloatType(0))) {
          if (eps_gt(a_ + a_ + e_ + e_ + f_, FloatType(0))) return false;
        }
        return true;
      }

    protected:

      // ---------------------------------------------------------------
      // Epsilon-aware comparison helpers
      // ---------------------------------------------------------------

      bool eps_lt(FloatType x, FloatType y) const
      {
        return x < y - epsilon_;
      }

      bool eps_gt(FloatType x, FloatType y) const
      {
        return y < x - epsilon_;
      }

      bool eps_eq(FloatType x, FloatType y) const
      {
        return !(eps_lt(x, y) || eps_lt(y, x));
      }

      // ---------------------------------------------------------------
      // Gruber sign analysis (epsilon-aware, matches reduction_base.def_test)
      // ---------------------------------------------------------------

      // Returns (n_zero, n_positive) for d_, e_, f_ using epsilon.
      af::tiny<int, 2>
      def_test() const
      {
        int n_zero = 0;
        int n_positive = 0;
        // eps_lt(0, x) means x > epsilon (positive in epsilon sense)
        // not eps_lt(x, 0) means x >= -epsilon (non-negative in epsilon sense)
        if (eps_lt(FloatType(0), d_)) n_positive++;
        else if (!eps_lt(d_, FloatType(0))) n_zero++;
        if (eps_lt(FloatType(0), e_)) n_positive++;
        else if (!eps_lt(e_, FloatType(0))) n_zero++;
        if (eps_lt(FloatType(0), f_)) n_positive++;
        else if (!eps_lt(f_, FloatType(0))) n_zero++;
        return af::tiny<int, 2>(n_zero, n_positive);
      }

      // True when all three are positive, or none are zero and exactly one positive.
      bool
      def_gt_0() const
      {
        af::tiny<int, 2> nz_np = def_test();
        return nz_np[1] == 3 || (nz_np[0] == 0 && nz_np[1] == 1);
      }

      // ---------------------------------------------------------------
      // Change-of-basis tracking
      // ---------------------------------------------------------------

      void
      cb_update(scitbx::mat3<IntFromFloatType> const& m)
      {
        if (n_iterations_ == iteration_limit_) {
          throw niggli_reduction_iteration_limit_exceeded(iteration_limit_);
        }
        r_inv_ = r_inv_ * m;
        n_iterations_ += 1;
      }

      // ---------------------------------------------------------------
      // N-actions (shared primitives; n3_true/false use epsilon)
      // ---------------------------------------------------------------

      void n1_action()
      {
        cb_update(scitbx::mat3<IntFromFloatType>(0,-1,0, -1,0,0, 0,0,-1));
        std::swap(a_, b_);
        std::swap(d_, e_);
      }

      void n2_action()
      {
        cb_update(scitbx::mat3<IntFromFloatType>(-1,0,0, 0,0,-1, 0,-1,0));
        std::swap(b_, c_);
        std::swap(e_, f_);
      }

      // Make d_,e_,f_ all non-negative (epsilon-aware sign check).
      void n3_true_action()
      {
        scitbx::mat3<IntFromFloatType> m(1);  // identity
        if (eps_lt(d_, FloatType(0))) m[0] = -1;
        if (eps_lt(e_, FloatType(0))) m[4] = -1;
        if (eps_lt(f_, FloatType(0))) m[8] = -1;
        cb_update(m);
        d_ = std::abs(d_);
        e_ = std::abs(e_);
        f_ = std::abs(f_);
      }

      // Make d_,e_,f_ all non-positive (epsilon-aware sign check).
      // z tracks which component is epsilon-zero (component index 0/1/2).
      void n3_false_action()
      {
        int f[3] = {1, 1, 1};
        int z = -1;
        if (eps_gt(d_, FloatType(0)))       f[0] = -1;
        else if (!eps_lt(d_, FloatType(0))) z = 0;
        if (eps_gt(e_, FloatType(0)))       f[1] = -1;
        else if (!eps_lt(e_, FloatType(0))) z = 1;
        if (eps_gt(f_, FloatType(0)))       f[2] = -1;
        else if (!eps_lt(f_, FloatType(0))) z = 2;
        if (f[0] * f[1] * f[2] < 0) {
          CCTBX_ASSERT(z != -1);
          f[z] = -1;
        }
        cb_update(scitbx::mat3<IntFromFloatType>(
          f[0],0,0, 0,f[1],0, 0,0,f[2]));
        d_ = -std::abs(d_);
        e_ = -std::abs(e_);
        f_ = -std::abs(f_);
      }

      // ---------------------------------------------------------------
      // A-actions (direct port of krivy_gruber_1976.py a*_action methods)
      // Note: a5/a6/a7 use exact sign comparison (d_>0, etc.), not epsilon.
      // This matches the Python source exactly.
      // ---------------------------------------------------------------

      void a1_action() { n1_action(); }
      void a2_action() { n2_action(); }
      void a3_action() { n3_true_action(); }
      void a4_action() { n3_false_action(); }

      void a5_action()
      {
        if (d_ > 0) {
          cb_update(scitbx::mat3<IntFromFloatType>(1,0,0, 0,1,-1, 0,0,1));
          c_ += b_ - d_;
          d_ -= b_ + b_;
          e_ -= f_;
        } else {
          cb_update(scitbx::mat3<IntFromFloatType>(1,0,0, 0,1,1, 0,0,1));
          c_ += b_ + d_;
          d_ += b_ + b_;
          e_ += f_;
        }
        CCTBX_ASSERT(c_ > 0);
      }

      void a6_action()
      {
        if (e_ > 0) {
          cb_update(scitbx::mat3<IntFromFloatType>(1,0,-1, 0,1,0, 0,0,1));
          c_ += a_ - e_;
          d_ -= f_;
          e_ -= a_ + a_;
        } else {
          cb_update(scitbx::mat3<IntFromFloatType>(1,0,1, 0,1,0, 0,0,1));
          c_ += a_ + e_;
          d_ += f_;
          e_ += a_ + a_;
        }
        CCTBX_ASSERT(c_ > 0);
      }

      void a7_action()
      {
        if (f_ > 0) {
          cb_update(scitbx::mat3<IntFromFloatType>(1,-1,0, 0,1,0, 0,0,1));
          b_ += a_ - f_;
          d_ -= e_;
          f_ -= a_ + a_;
        } else {
          cb_update(scitbx::mat3<IntFromFloatType>(1,1,0, 0,1,0, 0,0,1));
          b_ += a_ + f_;
          d_ += e_;
          f_ += a_ + a_;
        }
        CCTBX_ASSERT(b_ > 0);
      }

      void a8_action()
      {
        cb_update(scitbx::mat3<IntFromFloatType>(1,0,1, 0,1,1, 0,0,1));
        c_ += a_ + b_ + d_ + e_ + f_;
        d_ += b_ + b_ + f_;
        e_ += a_ + a_ + f_;
        CCTBX_ASSERT(c_ > 0);
      }

      // ---------------------------------------------------------------
      // Main step function — exact port of krivy_gruber_1976.step()
      // ---------------------------------------------------------------

      bool step()
      {
        // A1
        if (eps_gt(a_, b_)
            || (eps_eq(a_, b_) && eps_gt(std::abs(d_), std::abs(e_)))) {
          a1_action();
        }
        // A2
        if (eps_gt(b_, c_)
            || (eps_eq(b_, c_) && eps_gt(std::abs(e_), std::abs(f_)))) {
          a2_action();
          return true;
        }
        // A3 / A4
        if (def_gt_0()) {
          a3_action();
        } else {
          a4_action();
        }
        // A5
        if (eps_gt(std::abs(d_), b_)
            || (eps_eq(d_, b_)  && eps_lt(e_ + e_, f_))
            || (eps_eq(d_, -b_) && eps_lt(f_, FloatType(0)))) {
          a5_action();
          return true;
        }
        // A6
        if (eps_gt(std::abs(e_), a_)
            || (eps_eq(e_, a_)  && eps_lt(d_ + d_, f_))
            || (eps_eq(e_, -a_) && eps_lt(f_, FloatType(0)))) {
          a6_action();
          return true;
        }
        // A7
        if (eps_gt(std::abs(f_), a_)
            || (eps_eq(f_, a_)  && eps_lt(d_ + d_, e_))
            || (eps_eq(f_, -a_) && eps_lt(e_, FloatType(0)))) {
          a7_action();
          return true;
        }
        // A8
        if (eps_lt(d_ + e_ + f_ + a_ + b_, FloatType(0))
            || (eps_eq(d_ + e_ + f_ + a_ + b_, FloatType(0))
                && eps_gt(a_ + a_ + e_ + e_ + f_, FloatType(0)))) {
          a8_action();
          return true;
        }
        return false;
      }

      // ---------------------------------------------------------------
      // Data members
      // ---------------------------------------------------------------

      std::size_t iteration_limit_;
      scitbx::mat3<IntFromFloatType> r_inv_;
      std::size_t n_iterations_;
      FloatType epsilon_;
      FloatType a_, b_, c_, d_, e_, f_;
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_NIGGLI_REDUCTION_H
