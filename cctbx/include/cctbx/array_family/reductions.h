// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_REDUCTIONS_H
#define CCTBX_ARRAY_FAMILY_REDUCTIONS_H

#include <cctbx/array_family/ref.h>

namespace cctbx { namespace af {

  template <typename ElementType1, typename AccessorType1,
            typename ElementType2, typename AccessorType2>
  int
  cmp(
    const_ref<ElementType1, AccessorType1> const& a1,
    const_ref<ElementType2, AccessorType2> const& a2)
  {
    std::size_t sz_min = (a1.size() < a2.size() ? a1.size() : a2.size());
    for(std::size_t i=0;i<sz_min;i++) {
      if (a1[i] < a2[i]) return -1;
      if (a2[i] < a1[i]) return 1;
    }
    if (a1.size() < a2.size()) return -1;
    if (a2.size() < a1.size()) return 1;
    return 0;
  }

  template <typename ElementType, typename AccessorType>
  int
  cmp(
    const_ref<ElementType, AccessorType> const& a1,
    ElementType const& a2)
  {
    for(std::size_t i=0;i<a1.size();i++) {
      if (a1[i] < a2) return -1;
      if (a2 < a1[i]) return 1;
    }
    return 0;
  }

  template <typename ElementType, typename AccessorType>
  int
  cmp(
    ElementType const& a1,
    const_ref<ElementType, AccessorType> const& a2)
  {
    for(std::size_t i=0;i<a2.size();i++) {
      if (a1 < a2[i]) return -1;
      if (a2[i] < a1) return 1;
    }
    return 0;
  }

  template <typename ElementType, typename AccessorType>
  std::size_t
  max_index(const_ref<ElementType, AccessorType> const& a)
  {
    std::size_t result;
    if (a.size() > 0) {
      result = 0;
      for(std::size_t i=1;i<a.size();i++) {
        if (a[result] < a[i]) result = i;
      }
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  std::size_t
  min_index(const_ref<ElementType, AccessorType> const& a)
  {
    std::size_t result;
    if (a.size() > 0) {
      result = 0;
      for(std::size_t i=1;i<a.size();i++) {
        if (a[result] > a[i]) result = i;
      }
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  max(const_ref<ElementType, AccessorType> const& a)
  {
    ElementType result;
    if (a.size() > 0) {
      result = a[0];
      for(std::size_t i=1;i<a.size();i++) {
        if (result < a[i]) result = a[i];
      }
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  min(const_ref<ElementType, AccessorType> const& a)
  {
    ElementType result;
    if (a.size() > 0) {
      result = a[0];
      for(std::size_t i=1;i<a.size();i++) {
        if (result > a[i]) result = a[i];
      }
    }
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  sum(const_ref<ElementType, AccessorType> const& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  sum_sq(const_ref<ElementType, AccessorType> const& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += a[i] * a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  product(const_ref<ElementType, AccessorType> const& a)
  {
    std::size_t sz = a.size();
    if (sz == 0) return 0;
    ElementType result = 1;
    for(std::size_t i=0;i<sz;i++) result *= a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean(const_ref<ElementType, AccessorType> const& a)
  {
    return sum(a) / a.size();
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean_sq(const_ref<ElementType, AccessorType> const& a)
  {
    return sum_sq(a) / a.size();
  }

  template <typename ElementTypeValues, typename AccessorTypeValues,
            typename ElementTypeWeights, typename AccessorTypeWeights>
  ElementTypeValues
  mean_weighted(
    const_ref<ElementTypeValues, AccessorTypeValues> const& values,
    const_ref<ElementTypeWeights, AccessorTypeWeights> const& weights)
  {
    if (values.size() != weights.size()) throw_range_error();
    if (!values.size()) return 0;
    ElementTypeValues sum_vw = 0;
    ElementTypeWeights sum_w = 0;
    for(std::size_t i=0;i<values.size();i++) {
      sum_vw += values[i] * weights[i];
      sum_w += weights[i];
    }
    return sum_vw / sum_w;
  }

  template <typename ElementTypeValues, typename AccessorTypeValues,
            typename ElementTypeWeights, typename AccessorTypeWeights>
  ElementTypeValues
  mean_sq_weighted(
    const_ref<ElementTypeValues, AccessorTypeValues> const& values,
    const_ref<ElementTypeWeights, AccessorTypeWeights> const& weights)
  {
    if (values.size() != weights.size()) throw_range_error();
    if (!values.size()) return 0;
    ElementTypeValues sum_vvw = 0;
    ElementTypeWeights sum_w = 0;
    for(std::size_t i=0;i<values.size();i++) {
      sum_vvw += values[i] * values[i] * weights[i];
      sum_w += weights[i];
    }
    return sum_vvw / sum_w;
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_REDUCTIONS_H
