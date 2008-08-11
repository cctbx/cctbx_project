#ifndef CCTBX_UCTBX_FAST_MINIMUM_REDUCTION_H
#define CCTBX_UCTBX_FAST_MINIMUM_REDUCTION_H

#include <cctbx/uctbx.h>
#include <cctbx/error.h>

#if !defined(CCTBX_UCTBX_FAST_MINIMUM_REDUCTION_OVERZEALOUS_OPTIMIZER)
# if defined(__GNUC__) && __GNUC__ == 2 && __GNUC_MINOR__ == 96
#  define CCTBX_UCTBX_FAST_MINIMUM_REDUCTION_OVERZEALOUS_OPTIMIZER
#  include <scitbx/serialization/base_256.h>
# endif
#endif

namespace cctbx { namespace uctbx {

  //! Specific exception to indicate failure of reduction algorithm.
  class error_iteration_limit_exceeded : public error
  {
    public:
      error_iteration_limit_exceeded()
      :
        error("Iteration limit exceeded.")
      {}
  };

  //! Specific exception to indicate failure of reduction algorithm.
  class error_degenerate_unit_cell_parameters : public error
  {
    public:
      error_degenerate_unit_cell_parameters()
      :
        error("Degenerate unit cell parameters.")
      {}
  };

  //! Helper function to disable optimization.
  float spoil_optimization(float);
  //! Helper function to disable optimization.
  double spoil_optimization(double);

  //! Fast minimum-lengths cell reduction.
  /*! Based on the algorithm of Gruber (1973), Acta Cryst. A29, 433-440.
      Tests for equality are removed in order to make the algorithm
      numerically more stable. However, in some cases the algorithm
      still oscillates. This situation is analyzed in the implementation
      below and the cycle is terminated after the action in the false
      branch of N3.

      Attention:
        If this class is instantiated with a FloatType other than
        float or double a corresponding spoil_optimization() helper
        function must be defined!
   */
  template <typename FloatType=double,
            typename IntFromFloatType=int>
  class fast_minimum_reduction
  {
    public:
      //! Default contructor. Some data members are not initialized!
      fast_minimum_reduction() {}

      //! Executes the reduction algorithm.
      fast_minimum_reduction(
        uctbx::unit_cell const& unit_cell,
        std::size_t iteration_limit=100,
        FloatType multiplier_significant_change_test=16,
        std::size_t min_n_no_significant_change=2)
      :
        multiplier_significant_change_test_(
          multiplier_significant_change_test),
        min_n_no_significant_change_(min_n_no_significant_change),
        iteration_limit_(iteration_limit),
        r_inv_(scitbx::mat3<IntFromFloatType>(1)),
        n_iterations_(0),
        n_no_significant_change_(0),
        termination_due_to_significant_change_test_(false)
      {
        uc_sym_mat3 const& g = unit_cell.metrical_matrix();
        a_ = g[0];
        b_ = g[1];
        c_ = g[2];
        d_ = 2*g[5];
        e_ = 2*g[4];
        f_ = 2*g[3];
        last_abc_significant_change_test_ = scitbx::vec3<FloatType>(
          -a_,-b_,-c_);
        while (step());
      }

      //! Parameterization of Gruber (1973).
      af::double6
      as_gruber_matrix() const
      {
        return af::double6(a_, b_, c_, d_, e_, f_);
      }

      //! Parameterization of Niggli.
      af::double6
      as_niggli_matrix() const
      {
        return af::double6(a_, b_, c_, d_/2, e_/2, f_/2);
      }

      //! Parameterization as used by uctbx::unit_cell::metrical_matrix().
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

      //! Value as passed to the constructor.
      std::size_t
      iteration_limit() const { return iteration_limit_; }

      //! Value as passed to the constructor.
      FloatType
      multiplier_significant_change_test() const
      {
        return multiplier_significant_change_test_;
      }

      //! Value as passed to the constructor.
      std::size_t
      min_n_no_significant_change() const
      {
        return min_n_no_significant_change_;
      }

      //! Change-of-basis matrix.
      /*! Compatible with uctbx::unit_cell::change_basis()
       */
      scitbx::mat3<IntFromFloatType> const&
      r_inv() const { return r_inv_; }

      //! Number of iterations (actions executed).
      std::size_t
      n_iterations() const { return n_iterations_; }

      /*! \brief Indicates if the reduction was terminated because
          it was detected that the unit cell lengths did no longer
          change significantly.
       */
      bool
      termination_due_to_significant_change_test() const
      {
        return termination_due_to_significant_change_test_;
      }

      //! Cell type according to International Tables for Crystallography.
      /*! Possible values: 1 or 2. The value 0 is not possible if
          the algorithm terminates.
       */
      int
      type() const
      {
        int n_positive = def_test()[1];
        if (n_positive == 3) return 1;
        if (n_positive == 0) return 2;
        return 0;
      }

    protected:
      // greatest integer which is not greater than x
      inline
      IntFromFloatType
      entier(FloatType const& x)
      {
        IntFromFloatType result = static_cast<IntFromFloatType>(x);
        if (x-result < 0) result--;
        if (!(x-result < 1)) result++; // work around rounding errors
        return result;
      }

      af::tiny<int, 2>
      def_test() const
      {
        int n_zero = 0;
        int n_positive = 0;
        if (0 < d_) n_positive += 1;
        else if (!(d_ < 0)) n_zero += 1;
        if (0 < e_) n_positive += 1;
        else if (!(e_ < 0)) n_zero += 1;
        if (0 < f_) n_positive += 1;
        else if (!(f_ < 0)) n_zero += 1;
        return af::tiny<int, 2>(n_zero, n_positive);
      }

      bool
      def_gt_0() const
      {
        af::tiny<int, 2> nz_np = def_test();
        return nz_np[1] == 3 || (nz_np[0] == 0 && nz_np[1] == 1);
      }

      bool
      step()
      {
        // N1
        if (b_ < a_) {
          n1_action();
        }
        // N2
        if (c_ < b_) {
          n2_action();
          return true;
        }
        // N3
        if (def_gt_0()) {
          n3_true_action();
        }
        else {
          n3_false_action();
          if (!significant_change_test()) {
            return false;
          }
        }
        if (b2_action()) return true;
        if (b3_action()) return true;
        if (b4_action()) return true;
        if (b5_action()) return true;
        return false;
      }

      FloatType
      significant_change_test(FloatType const& new_value, std::size_t i) const
      {
        FloatType m_new = multiplier_significant_change_test_ * new_value;
        FloatType diff = new_value - last_abc_significant_change_test_[i];
        FloatType m_new_plus_diff = spoil_optimization(m_new + diff);
#if defined(CCTBX_UCTBX_FAST_MINIMUM_REDUCTION_OVERZEALOUS_OPTIMIZER)
        {
          char buf[4*sizeof(FloatType)];
          scitbx::serialization::base_256::floating_point
            ::to_string(buf, m_new_plus_diff);
          scitbx::serialization::base_256::floating_point
            ::from_string<FloatType> fs(buf);
          m_new_plus_diff = fs.value;
        }
#endif
        FloatType m_new_plus_diff_minus_m_new = m_new_plus_diff - m_new;
        return m_new_plus_diff_minus_m_new != 0;
      }

      bool
      significant_change_test()
      {
        if (   significant_change_test(a_, 0)
            || significant_change_test(b_, 1)
            || significant_change_test(c_, 2)) {
          n_no_significant_change_ = 0;
        }
        else {
          n_no_significant_change_++;
          if (n_no_significant_change_ == min_n_no_significant_change_) {
            return false;
          }
        }
        last_abc_significant_change_test_ = scitbx::vec3<FloatType>(a_,b_,c_);
        return true;
      }

      void
      cb_update(scitbx::mat3<IntFromFloatType> const& m)
      {
        if (n_iterations_ == iteration_limit_) {
          throw error_iteration_limit_exceeded();
        }
        r_inv_ = r_inv_ * m;
        n_iterations_ += 1;
      }

      void
      n1_action()
      {
        cb_update(scitbx::mat3<IntFromFloatType>(0,-1,0, -1,0,0, 0,0,-1));
        std::swap(a_, b_);
        std::swap(d_, e_);
      }

      void
      n2_action()
      {
        cb_update(scitbx::mat3<IntFromFloatType>(-1,0,0, 0,0,-1, 0,-1,0));
        std::swap(b_, c_);
        std::swap(e_, f_);
      }

      void
      n3_true_action()
      {
        scitbx::mat3<IntFromFloatType> m(1);
        if (d_ < 0) m[0] = -1;
        if (e_ < 0) m[4] = -1;
        if (f_ < 0) m[8] = -1;
        cb_update(m);
        d_ = scitbx::fn::absolute(d_);
        e_ = scitbx::fn::absolute(e_);
        f_ = scitbx::fn::absolute(f_);
      }

      void
      n3_false_action()
      {
        scitbx::mat3<IntFromFloatType> m(1);
        int z = -1;
        if (0 < d_) m[0] = -1;
        else if (!(d_ < 0)) z = 0;
        if (0 < e_) m[4] = -1;
        else if (!(e_ < 0)) z = 4;
        if (0 < f_) m[8] = -1;
        else if (!(f_ < 0)) z = 8;
        if (m[0]*m[4]*m[8] < 0) {
          CCTBX_ASSERT(z != -1);
          m[z] = -1;
        }
        cb_update(m);
        d_ = -scitbx::fn::absolute(d_);
        e_ = -scitbx::fn::absolute(e_);
        f_ = -scitbx::fn::absolute(f_);
      }

      bool
      b2_action()
      {
        if (!(b_ < scitbx::fn::absolute(d_))) return false;
        IntFromFloatType j = entier((d_+b_)/(2*b_));
        if (j == 0) return false;
        cb_update(scitbx::mat3<IntFromFloatType>(1,0,0,0,1,-j,0,0,1));
        c_ += j*j*b_ - j*d_;
        d_ -= 2*j*b_;
        e_ -= j*f_;
        if (!(0 < c_)) throw error_degenerate_unit_cell_parameters();
        return true;
      }

      bool
      b3_action()
      {
        if (!(a_ < scitbx::fn::absolute(e_))) return false;
        IntFromFloatType j = entier((e_+a_)/(2*a_));
        if (j == 0) return false;
        cb_update(scitbx::mat3<IntFromFloatType>(1,0,-j,0,1,0,0,0,1));
        c_ += j*j*a_ - j*e_;
        d_ -= j*f_;
        e_ -= 2*j*a_;
        if (!(0 < c_)) throw error_degenerate_unit_cell_parameters();
        return true;
      }

      bool
      b4_action()
      {
        if (!(a_ < scitbx::fn::absolute(f_))) return false;
        IntFromFloatType j = entier((f_+a_)/(2*a_));
        if (j == 0) return false;
        cb_update(scitbx::mat3<IntFromFloatType>(1,-j,0,0,1,0,0,0,1));
        b_ += j*j*a_ - j*f_;
        d_ -= j*e_;
        f_ -= 2*j*a_;
        if (!(0 < b_)) throw error_degenerate_unit_cell_parameters();
        return true;
      }

      bool
      b5_action()
      {
        FloatType de = d_ + e_;
        FloatType fab = f_ + a_ + b_;
        if (!(de+fab < 0)) return false;
        IntFromFloatType j = entier((de+fab)/(2*fab));
        if (j == 0) return false;
        cb_update(scitbx::mat3<IntFromFloatType>(1,0,-j,0,1,-j,0,0,1));
        c_ += j*j*fab-j*de;
        d_ -= j*(2*b_+f_);
        e_ -= j*(2*a_+f_);
        if (!(0 < c_)) throw error_degenerate_unit_cell_parameters();
        return true;
      }

      FloatType multiplier_significant_change_test_;
      std::size_t min_n_no_significant_change_;
      std::size_t iteration_limit_;
      scitbx::mat3<IntFromFloatType> r_inv_;
      std::size_t n_iterations_;
      std::size_t n_no_significant_change_;
      bool termination_due_to_significant_change_test_;
      FloatType a_, b_, c_, d_, e_, f_;
      scitbx::vec3<FloatType> last_abc_significant_change_test_;
  };

}} // namespace cctbx::uctbx

#endif // CCTBX_UCTBX_FAST_MINIMUM_REDUCTION_H
