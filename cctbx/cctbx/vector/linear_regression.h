// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_LINEAR_REGRESSION_H
#define CCTBX_VECTOR_LINEAR_REGRESSION_H

#include <cctbx/error.h>

namespace cctbx { namespace vector {

  template <class IntegerType, class FloatType>
  struct linear_regression_core
  {
    void reset() { b = m = cc = 0; }
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
        b = min_y;
        return;
      }
      FloatType fn(n);
      FloatType dx = std::max(f_abs(min_x - sum_x / fn),
                              f_abs(max_x - sum_x / fn));
      FloatType dy = std::max(f_abs(min_y - sum_y / fn),
                              f_abs(max_y - sum_y / fn));
      if (dx == 0) return;
      if (dy == 0) {
        b = sum_y / fn;
        return;
      }
      if (dx < dy * epsilon) return;
      if (dy < dx * epsilon) {
        b = sum_y / fn;
        return;
      }
      FloatType d = fn * sum_x2 - sum_x * sum_x;
      if (d != 0) {
        b = (sum_x2 * sum_y - sum_x * sum_xy) / d;
        m = (fn * sum_xy - sum_x * sum_y) / d;
      }
      d =   (sum_x2 - sum_x * sum_x / fn)
          * (sum_y2 - sum_y * sum_y / fn);
      if (d > 0.)
        cc = (sum_xy - sum_x * sum_y / fn) / std::sqrt(d);
    }
    //! Work-around for broken compilers.
    FloatType f_abs(const FloatType& x) {
      if (x < 0) return -x;
      return x;
    }
    //! Linear regression coefficients b: y = mx + b.
    FloatType b;
    //! Linear regression coefficients m: y = mx + b.
    FloatType m;
    //! Standard linear correlation coefficient.
    FloatType cc;
  };

  template <class VectorType>
  struct linear_regression
    : linear_regression_core<typename VectorType::size_type,
                             typename VectorType::value_type>
  {
    typedef typename VectorType::size_type size_type;
    typedef typename VectorType::value_type value_type;

    linear_regression() {}
    linear_regression(const VectorType& x,
                      const VectorType& y,
                      const value_type& epsilon = 1.e-6)
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
      set(n, min_x, max_x, min_y, max_y, sum_x, sum_x2, sum_y, sum_y2, sum_xy,
          epsilon);
    }
  };

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_LINEAR_REGRESSION_H
