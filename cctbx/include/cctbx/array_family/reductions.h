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

  template <typename ElementType, typename AccessorType>
  std::size_t
  max_index(const const_ref<ElementType, AccessorType>& a)
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
  min_index(const const_ref<ElementType, AccessorType>& a)
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
  max(const const_ref<ElementType, AccessorType>& a)
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
  min(const const_ref<ElementType, AccessorType>& a)
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
  sum(const const_ref<ElementType, AccessorType>& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  sum_sq(const const_ref<ElementType, AccessorType>& a)
  {
    ElementType result = 0;
    for(std::size_t i=0;i<a.size();i++) result += a[i] * a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  product(const const_ref<ElementType, AccessorType>& a)
  {
    ElementType result = 1;
    for(std::size_t i=0;i<a.size();i++) result *= a[i];
    return result;
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean(const const_ref<ElementType, AccessorType>& a)
  {
    return sum(a) / a.size();
  }

  template <typename ElementType, typename AccessorType>
  ElementType
  mean_sq(const const_ref<ElementType, AccessorType>& a)
  {
    return sum_sq(a) / a.size();
  }

  template <typename ElementTypeValues, typename AccessorTypeValues,
            typename ElementTypeWeights, typename AccessorTypeWeights>
  ElementTypeValues
  mean_weighted(
    const const_ref<ElementTypeValues, AccessorTypeValues>& values,
    const const_ref<ElementTypeWeights, AccessorTypeWeights>& weights)
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
    const const_ref<ElementTypeValues, AccessorTypeValues>& values,
    const const_ref<ElementTypeWeights, AccessorTypeWeights>& weights)
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
