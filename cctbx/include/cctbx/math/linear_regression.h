// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: modified version of cctbx/vector/linear_regresion.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MATH_LINEAR_REGRESSION_H
#define CCTBX_MATH_LINEAR_REGRESSION_H

#include <cctbx/error.h>
#include <cctbx/fixes/cmath>
#include <cctbx/array_family/tiny.h>

namespace cctbx { namespace math {

  template <class IntegerType, class FloatType>
  class linear_regression_core
  {
    public:
      void reset() {
        m_is_well_defined = false;
        m_b = m_m = m_cc = FloatType(0);
      }
      void set(const IntegerType& n,
               const FloatType& min_x, const FloatType& max_x,
               const FloatType& min_y, const FloatType& max_y,
               const FloatType& sum_x, const FloatType& sum_x2,
               const FloatType& sum_y, const FloatType& sum_y2,
               const FloatType& sum_xy,
               const FloatType& epsilon = 1.e-6)
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
        d =   (sum_x2 - sum_x * sum_x / fn)
            * (sum_y2 - sum_y * sum_y / fn);
        if (d > 0.)
          m_cc = (sum_xy - sum_x * sum_y / fn) / std::sqrt(d);
        m_is_well_defined = true;
      }
      //! Work-around for broken compilers.
      static FloatType f_abs(const FloatType& x) {
        if (x < 0) return -x;
        return x;
      }
      //! Flag.
      bool is_well_defined() const { return m_is_well_defined; }
      //! Linear regression coefficients b: y = mx + b.
      const FloatType& b() const { return m_b; }
      //! Linear regression coefficients m: y = mx + b.
      const FloatType& m() const { return m_m; }
      //! Standard linear correlation coefficient.
      const FloatType& cc() const { return m_cc; }
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
      linear_regression(const af::const_ref<FloatType>& x,
                        const af::const_ref<FloatType>& y,
                        const FloatType& epsilon = 1.e-6)
      {
        size_type n = x.size();
        cctbx_assert(n == y.size());
        if (n == 0) {
          reset();
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
        set(n, min_x, max_x, min_y, max_y,
            sum_x, sum_x2, sum_y, sum_y2, sum_xy,
            epsilon);
      }
  };

}} // namespace cctbx::math

#endif // CCTBX_MATH_LINEAR_REGRESSION_H
