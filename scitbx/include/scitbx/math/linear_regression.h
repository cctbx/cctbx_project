/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copy of cctbx/math/linear_regression.h (rwgk)
     2002 Mar: modified version of cctbx/vector/linear_regresion.h (rwgk)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_MATH_LINEAR_REGRESSION_H
#define SCITBX_MATH_LINEAR_REGRESSION_H

#include <cmath>
#include <scitbx/error.h>
#include <scitbx/array_family/tiny.h>

namespace scitbx { namespace math {

  template <typename FloatType>
  FloatType
  std_linear_correlation(
    FloatType const& sum_weights,
    FloatType const& sum_x, FloatType const& sum_x2,
    FloatType const& sum_y, FloatType const& sum_y2,
    FloatType const& sum_xy)
  {
    FloatType cc_denom2 =   (sum_x2 - sum_x * sum_x / sum_weights)
                          * (sum_y2 - sum_y * sum_y / sum_weights);
    FloatType cc(1);
    if (cc_denom2 > FloatType(0)) {
      cc = (sum_xy - sum_x * sum_y / sum_weights) / std::sqrt(cc_denom2);
    }
    return cc;
  }

  template <class IntegerType, class FloatType>
  class linear_regression_core
  {
    public:
      void reset() {
        m_is_well_defined = false;
        m_b = m_m = m_cc = FloatType(0);
      }
      void set(IntegerType const& n,
               FloatType const& min_x, FloatType const& max_x,
               FloatType const& min_y, FloatType const& max_y,
               FloatType const& sum_x, FloatType const& sum_x2,
               FloatType const& sum_y, FloatType const& sum_y2,
               FloatType const& sum_xy,
               FloatType const& epsilon = 1.e-6)
      {
        reset();
        if (n < 1) return;
        if (min_x == max_x) return;
        if (min_y == max_y) {
          m_b = min_y;
          m_cc = FloatType(1);
          m_is_well_defined = true;
          return;
        }
        FloatType fn(n);
        FloatType dx = std::max(f_abs(min_x - sum_x / fn),
                                f_abs(max_x - sum_x / fn));
        FloatType dy = std::max(f_abs(min_y - sum_y / fn),
                                f_abs(max_y - sum_y / fn));
        if (dx == 0) return;
        if (dy == 0) {
          m_b = sum_y / fn;
          m_cc = FloatType(1);
          m_is_well_defined = true;
          return;
        }
        if (dx < dy * epsilon) return;
        FloatType d = fn * sum_x2 - sum_x * sum_x;
        if (d != 0) {
          m_b = (sum_x2 * sum_y - sum_x * sum_xy) / d;
          m_m = (fn * sum_xy - sum_x * sum_y) / d;
        }
        m_cc = std_linear_correlation(
          fn, sum_x, sum_x2, sum_y, sum_y2, sum_xy);
        m_is_well_defined = true;
      }
      //! Work-around for broken compilers.
      static FloatType f_abs(FloatType const& x) {
        if (x < 0) return -x;
        return x;
      }
      //! Flag.
      bool is_well_defined() const { return m_is_well_defined; }
      //! Linear regression coefficients b: y = mx + b.
      FloatType const& b() const { return m_b; }
      //! Linear regression coefficients m: y = mx + b.
      FloatType const& m() const { return m_m; }
      //! Standard linear correlation coefficient.
      FloatType const& cc() const { return m_cc; }
    protected:
      bool m_is_well_defined;
      FloatType m_b;
      FloatType m_m;
      FloatType m_cc;
  };

  template <class FloatType>
  class linear_regression
    : public linear_regression_core<std::size_t, FloatType>
  {
    public:
      typedef std::size_t size_type;
      typedef FloatType value_type;

      linear_regression() {}
      linear_regression(af::const_ref<FloatType> const& x,
                        af::const_ref<FloatType> const& y,
                        FloatType const& epsilon = 1.e-6)
      {
        size_type n = x.size();
        SCITBX_ASSERT(n == y.size());
        if (n == 0) {
          this->reset();
          return;
        }
        value_type min_x = x[0];
        value_type max_x = x[0];
        value_type min_y = y[0];
        value_type max_y = y[0];
        value_type sum_x = x[0];
        value_type sum_x2 = x[0] * x[0];
        value_type sum_y = y[0];
        value_type sum_y2 = y[0] * y[0];
        value_type sum_xy = x[0] * y[0];
        for(size_type i=1;i<n;i++) {
          if (min_x > x[i]) min_x = x[i];
          if (max_x < x[i]) max_x = x[i];
          if (min_y > y[i]) min_y = y[i];
          if (max_y < y[i]) max_y = y[i];
          sum_x += x[i];
          sum_x2 += x[i] * x[i];
          sum_y += y[i];
          sum_y2 += y[i] * y[i];
          sum_xy += x[i] * y[i];
        }
        this->set(n, min_x, max_x, min_y, max_y,
                  sum_x, sum_x2, sum_y, sum_y2, sum_xy,
                  epsilon);
      }
  };

}} // namespace scitbx::math

#endif // SCITBX_MATH_LINEAR_REGRESSION_H
