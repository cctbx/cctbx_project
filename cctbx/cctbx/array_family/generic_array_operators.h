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
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a++, result++) {
        *result = op(*a);
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
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
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a++, result++) {
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

  // inplace_op(array, array)
  template <typename InPlaceBinaryOperatorType,
            typename ElementType1,
            typename ElementType2>
  inline void
  array_operation_in_place_binary(
    const InPlaceBinaryOperatorType& op,
    ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz)
  {
    ElementType1* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++, a2++) op(*a1, ElementType1(*a2));
  }

  // inplace_op(array, scalar)
  template <typename InPlaceBinaryOperatorType,
            typename ElementType>
  inline void
  array_operation_in_place_binary(
    const InPlaceBinaryOperatorType& op,
    ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz)
  {
    ElementType* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++) op(*a1, a2);
  }

  template <typename BooleanOperatorType,
            typename ElementType1,
            typename ElementType2>
  inline bool
  array_operation_reducing_boolean_op(
    const BooleanOperatorType& op,
    const ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz1,
    const std::size_t& sz2)
  {
    if (sz1 != sz2) throw_range_error();
    const ElementType1* a1_end = a1 + sz1;
    for(;a1 != a1_end; a1++, a2++) {
      if (!op(*a1, *a2)) return false;
    }
    return true;
  }

  template <typename BooleanOperatorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_op(
    const BooleanOperatorType& op,
    const ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz1)
  {
    const ElementType* a1_end = a1 + sz1;
    for(;a1 != a1_end; a1++) {
      if (!op(*a1, a2)) return false;
    }
    return true;
  }

  template <typename BooleanOperatorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_op(
    const BooleanOperatorType& op,
    const ElementType& a1,
    const ElementType* a2,
    const std::size_t& sz2)
  {
    const ElementType* a2_end = a2 + sz2;
    for(;a2 != a2_end; a2++) {
      if (!op(a1, *a2)) return false;
    }
    return true;
  }

  template <typename BooleanOperatorType,
            typename ElementType1,
            typename ElementType2>
  inline bool
  array_operation_reducing_boolean_op_not_equal_to(
    const BooleanOperatorType& op,
    const ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz1,
    const std::size_t& sz2)
  {
    if (sz1 != sz2) throw_range_error();
    const ElementType1* a1_end = a1 + sz1;
    for(;a1 != a1_end; a1++, a2++) {
      if (op(*a1, *a2)) return true;
    }
    return false;
  }

  template <typename BooleanOperatorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_op_not_equal_to(
    const BooleanOperatorType& op,
    const ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz1)
  {
    const ElementType* a1_end = a1 + sz1;
    for(;a1 != a1_end; a1++) {
      if (op(*a1, a2)) return true;
    }
    return false;
  }

  template <typename BooleanOperatorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_op_not_equal_to(
    const BooleanOperatorType& op,
    const ElementType& a1,
    const ElementType* a2,
    const std::size_t& sz2)
  {
    const ElementType* a2_end = a2 + sz2;
    for(;a2 != a2_end; a2++) {
      if (op(a1, *a2)) return true;
    }
    return false;
  }

  template <typename BooleanOperatorType,
            typename ElementType1,
            typename ElementType2>
  inline bool
  array_operation_reducing_boolean_op_greater_less(
    const BooleanOperatorType& op,
    const ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz1,
    const std::size_t& sz2)
  {
    if (sz1 != sz2) throw_range_error();
    const ElementType1* a1_end = a1 + sz1;
    for(;a1 != a1_end; a1++, a2++) {
      if (op(*a1, *a2)) return true;
      if (op(*a2, *a1)) return false;
    }
    return false;
  }

  template <typename BooleanOperatorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_op_greater_less(
    const BooleanOperatorType& op,
    const ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz1)
  {
    const ElementType* a1_end = a1 + sz1;
    for(;a1 != a1_end; a1++) {
      if (op(*a1, a2)) return true;
      if (op(a2, *a1)) return false;
    }
    return false;
  }

  template <typename BooleanOperatorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_op_greater_less(
    const BooleanOperatorType& op,
    const ElementType& a1,
    const ElementType* a2,
    const std::size_t& sz2)
  {
    const ElementType* a2_end = a2 + sz2;
    for(;a2 != a2_end; a2++) {
      if (op(a1, *a2)) return true;
      if (op(*a2, a1)) return false;
    }
    return false;
  }

#define CCTBX_ARRAY_FAMILY_ARRAY_OPERATION_1ARG_ELEMENT_WISE(f) \
  { \
    ElementType* result_start = result.begin(); \
    ElementType* result_iter = result_start; \
    try { \
      ElementType* result_end = result_iter + a.size(); \
      const ElementType* a_iter = a.begin(); \
      for(;result_iter != result_end; a_iter++, result_iter++) { \
        *result_iter = f(*a_iter); \
      } \
    } \
    catch (...) { \
      cctbx::af::detail::destroy_array_elements(result_start, result_iter, \
        has_trivial_destructor<ElementType>::value()); \
      throw; \
    } \
  }

#define CCTBX_ARRAY_FAMILY_ARRAY_OPERATION_2ARG_ELEMENT_WISE_A_A(f) \
  { \
    ElementType1* result_start = result.begin(); \
    ElementType1* result_iter = result_start; \
    try { \
      ElementType1* result_end = result_iter + a1.size(); \
      const ElementType1* a1_iter = a1.begin(); \
      const ElementType2* a2_iter = a2.begin(); \
      for(;result_iter != result_end; a1_iter++, a2_iter_++, result_iter++) { \
        *result_iter = f(*a1_iter, *a2_iter); \
      } \
    } \
    catch (...) { \
      cctbx::af::detail::destroy_array_elements(result_start, result_iter, \
        has_trivial_destructor<ElementType>::value()); \
      throw; \
    } \
  }

#define CCTBX_ARRAY_FAMILY_ARRAY_OPERATION_2ARG_ELEMENT_WISE_A_S(f) \
  { \
    ElementType* result_start = result.begin(); \
    ElementType* result_iter = result_start; \
    try { \
      ElementType* result_end = result_iter + a1.size(); \
      const ElementType* a1_iter = a1.begin(); \
      for(;result_iter != result_end; a1_iter++, result_iter++) { \
        *result_iter = f(*a1_iter, a2); \
      } \
    } \
    catch (...) { \
      cctbx::af::detail::destroy_array_elements(result_start, result_iter, \
        has_trivial_destructor<ElementType>::value()); \
      throw; \
    } \
  }

#define CCTBX_ARRAY_FAMILY_ARRAY_OPERATION_2ARG_ELEMENT_WISE_S_A(f) \
  { \
    ElementType* result_start = result.begin(); \
    ElementType* result_iter = result_start; \
    try { \
      ElementType* result_end = result_iter + a2.size(); \
      const ElementType* a2_iter = a2.begin(); \
      for(;result_iter != result_end; a2_iter++, result_iter++) { \
        *result_iter = f(a1, *a2_iter); \
      } \
    } \
    catch (...) { \
      cctbx::af::detail::destroy_array_elements(result_start, result_iter, \
        has_trivial_destructor<ElementType>::value()); \
      throw; \
    } \
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H
