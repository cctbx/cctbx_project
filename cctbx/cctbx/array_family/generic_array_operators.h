// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H
#define CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H

#include <cctbx/array_family/misc.h>

namespace cctbx { namespace af {

  template <typename UnaryOperatorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_unary(
    const UnaryOperatorType& op,
    const ElementType* a,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      for(;a != a + sz; a++, result++) {
        *result = op(*a);
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result);
      throw;
    }
  }

  template <typename UnaryOperatorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_unary(
    const UnaryOperatorType& op,
    const ElementType* a,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    for(;a != a + sz; a++, result++) {
      *result = op(*a);
    }
  }

  // op(array, array), non-POD result type
  template <typename BinaryOperatorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary(
    const BinaryOperatorType& op,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, a2++, result++) {
        *result = op(ElementTypeResult(*a1), ElementTypeResult(*a2));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result);
      throw;
    }
  }

  // op(array, array), POD result type
  template <typename BinaryOperatorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_binary(
    const BinaryOperatorType& op,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, a2++, result++) {
      *result = op(ElementTypeResult(*a1), ElementTypeResult(*a2));
    }
  }

  // op(array, scalar), non-POD result type
  template <typename BinaryOperatorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary(
    const BinaryOperatorType& op,
    const ElementType* a1,
    const ElementType& a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, result++) {
        *result = op(ElementTypeResult(*a1), ElementTypeResult(a2));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result);
      throw;
    }
  }

  // op(array, scalar), POD result type
  template <typename BinaryOperatorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_binary(
    const BinaryOperatorType& op,
    const ElementType& a1,
    const ElementType* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, result++) {
      *result = op(ElementTypeResult(*a1), ElementTypeResult(a2));
    }
  }

  // op(scalar, array), non-POD result type
  template <typename BinaryOperatorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary(
    const BinaryOperatorType& op,
    const ElementType& a1,
    const ElementType* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a2++, result++) {
        *result = op(ElementTypeResult(a1), ElementTypeResult(*a2));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result);
      throw;
    }
  }

  // op(scalar, array), POD result type
  template <typename BinaryOperatorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_binary(
    const BinaryOperatorType& op,
    const ElementType* a1,
    const ElementType& a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a2++, result++) {
      *result = op(ElementTypeResult(a1), ElementTypeResult(*a2));
    }
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H
