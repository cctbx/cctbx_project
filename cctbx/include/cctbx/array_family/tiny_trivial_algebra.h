// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: modified copy of parts of matrixlite.h (rwgk)
     2001 Oct 16: Moved tensor transformations from adptbx (rwgk)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 May 07 added: identidy, isDiagonal, transpose
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILAY_TINY_TRIVIAL_ALGEBRA_H
#define CCTBX_ARRAY_FAMILAY_TINY_TRIVIAL_ALGEBRA_H

#include <cctbx/array_family/tiny.h>
#include <cctbx/array_family/misc_functions.h>

namespace cctbx { namespace af {

  template<typename AnyType, std::size_t N>
  inline
  bool
  operator==(const af::tiny<AnyType,N>& x, const af::tiny<AnyType,N>& y) {
      for (std::size_t i = 0; i < x.size(); i++)
          if (x[i] != y[i]) return false;
      return true;
  }

  template<typename AnyType, std::size_t N>
  inline
  bool
  operator==(const af::tiny<AnyType,N>& x, const AnyType& value) {
      for (std::size_t i = 0; i < x.size(); i++)
          if (x[i] != value) return false;
      return true;
  }

  template<typename AnyType, std::size_t N>
  inline
  bool
  operator!=(const af::tiny<AnyType,N>& x, const af::tiny<AnyType,N>& y) {
      for (std::size_t i = 0; i < x.size(); i++)
          if (x[i] != y[i]) return true;
      return false;
  }

  template<typename AnyType, std::size_t N>
  inline
  bool
  operator!=(const af::tiny<AnyType,N>& x, const AnyType& value) {
      for (std::size_t i = 0; i < x.size(); i++)
          if (x[i] != value) return true;
      return false;
  }

  template <typename AnyType, std::size_t N>
  inline
  std::size_t
  min_index(const af::tiny<AnyType, N>& a) {
    return std::min_element(a.begin(), a.end()) - a.begin();
  }

  template <typename AnyType, std::size_t N>
  inline
  std::size_t
  max_index(const af::tiny<AnyType, N>& a) {
    return std::max_element(a.begin(), a.end()) - a.begin();
  }

  template<typename NumType, std::size_t N>
  inline
  af::tiny<NumType,N>
  operator+(const af::tiny<NumType,N>& lhs,
            const af::tiny<NumType,N>& rhs) {
      af::tiny<NumType,N> result;
      for (std::size_t i = 0; i < lhs.size(); i++) {
          result[i] = lhs[i] + rhs[i];
      }
      return result;
  }

  template<typename NumType, std::size_t N>
  inline
  af::tiny<NumType,N>
  operator-(const af::tiny<NumType,N>& lhs,
            const af::tiny<NumType,N>& rhs) {
      af::tiny<NumType,N> result;
      for (std::size_t i = 0; i < lhs.size(); i++) {
          result[i] = lhs[i] - rhs[i];
      }
      return result;
  }

  template<typename NumType, std::size_t N>
  inline
  af::tiny<NumType,N>
  operator-(const af::tiny<NumType,N>& rhs) {
      af::tiny<NumType,N> result;
      for (std::size_t i = 0; i < rhs.size(); i++) {
          result[i] = -rhs[i];
      }
      return result;
  }

  template<typename NumType, std::size_t N>
  inline
  af::tiny<NumType,N>
  operator*(const af::tiny<NumType,N>& lhs,
            const af::tiny<NumType,N>& rhs) {
      af::tiny<NumType,N> result;
      for (std::size_t i = 0; i < lhs.size(); i++) {
          result[i] = lhs[i] * rhs[i];
      }
      return result;
  }

  template<typename NumType, std::size_t N>
  inline
  af::tiny<NumType,N>
  operator*(const NumType& lhs,
            const af::tiny<NumType,N>& rhs) {
      af::tiny<NumType,N> result;
      for (std::size_t i = 0; i < rhs.size(); i++) result[i] = lhs * rhs[i];
      return result;
  }

  template<typename NumType, std::size_t N>
  inline
  af::tiny<NumType,N>
  operator/(const af::tiny<NumType,N>& lhs,
            const NumType& rhs) {
      af::tiny<NumType,N> result;
      for (std::size_t i = 0; i < lhs.size(); i++) result[i] = lhs[i] / rhs;
      return result;
  }

  template<typename NumType, std::size_t N>
  inline
  NumType
  sum(const af::tiny<NumType,N>& a) {
    NumType result = 0;
    for (std::size_t i = 0; i < a.size(); i++) result += a[i];
    return result;
  }

  template <typename NumType, std::size_t N>
  inline
  af::tiny<NumType, N>
  abs(const af::tiny<NumType, N>& a)
  {
    af::tiny<NumType, N> result;
    for (std::size_t i = 0; i < N; i++) {
      if (a[i] < 0) result[i] = -a[i];
      else          result[i] =  a[i];
    }
    return result;
  }

  template <typename FloatType, std::size_t N>
  inline
  tiny<bool, N>
  approx_equal_scaled(const af::tiny<FloatType, N>& a,
                      const af::tiny<FloatType, N>& b,
                      FloatType scaled_tolerance) {
    tiny<bool, N> result;
    for (std::size_t i = 0; i < N; i++) {
      result[i] = cctbx::af::approx_equal_scaled(a[i], b[i], scaled_tolerance);
    }
    return result;
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILAY_TINY_TRIVIAL_ALGEBRA_H
