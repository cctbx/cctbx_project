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

  template <typename UnaryFunctorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_unary(
    const UnaryFunctorType& ftor,
    const ElementType* a,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a++, result++) {
        new (result) ElementTypeResult(ftor(*a));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  template <typename UnaryFunctorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_unary(
    const UnaryFunctorType& ftor,
    const ElementType* a,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a++, result++) {
      *result = ftor(*a);
    }
  }

  // ftor(array, array), non-POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary(
    const BinaryFunctorType& ftor,
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
        new (result) ElementTypeResult(ftor(*a1, *a2));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // ftor(array, array), POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_binary(
    const BinaryFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, a2++, result++) {
      *result = ftor(*a1, *a2);
    }
  }

  // ftor(array, scalar), non-POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary(
    const BinaryFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2& a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, result++) {
        new (result) ElementTypeResult(ftor(*a1, a2));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // ftor(array, scalar), POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_binary(
    const BinaryFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2& a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, result++) {
      *result = ftor(*a1, a2);
    }
  }

  // ftor(scalar, array), non-POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary(
    const BinaryFunctorType& ftor,
    const ElementType1& a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a2++, result++) {
        new (result) ElementTypeResult(ftor(a1, *a2));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // ftor(scalar, array), POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_binary(
    const BinaryFunctorType& ftor,
    const ElementType1& a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a2++, result++) {
      *result = ftor(a1, *a2);
    }
  }

  // functor(array, array, addl), non-POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary_addl(
    const BinaryFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    const ElementType3& a3,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, a2++, result++) {
        new (result) ElementTypeResult(ftor(*a1, *a2, a3));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array, array, addl), POD result type
  template <typename BinaryFunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  inline void
  array_operation_binary_addl(
    const BinaryFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    const ElementType3& a3,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, a2++, result++) {
      *result = ftor(*a1, *a2, a3);
    }
  }

  // functor(array, scalar), non-POD result type
  template <typename BinaryFunctorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary_addl(
    const BinaryFunctorType& ftor,
    const ElementType* a1,
    const ElementType& a2,
    const ElementType& a3,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, result++) {
        new (result) ElementTypeResult(ftor(*a1, a2, a3));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array, scalar), POD result type
  template <typename BinaryFunctorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_binary_addl(
    const BinaryFunctorType& ftor,
    const ElementType* a1,
    const ElementType& a2,
    const ElementType& a3,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, result++) {
      *result = ftor(*a1, a2, a3);
    }
  }

  // functor(scalar, array), non-POD result type
  template <typename BinaryFunctorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_binary_addl(
    const BinaryFunctorType& ftor,
    const ElementType& a1,
    const ElementType* a2,
    const ElementType& a3,
    ElementTypeResult* result,
    const std::size_t& sz,
    false_type) // non-POD result type
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a2++, result++) {
        new (result) ElementTypeResult(ftor(a1, *a2, a3));
      }
    }
    catch (...) {
      cctbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(scalar, array), POD result type
  template <typename BinaryFunctorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_binary_addl(
    const BinaryFunctorType& ftor,
    const ElementType& a1,
    const ElementType* a2,
    const ElementType& a3,
    ElementTypeResult* result,
    const std::size_t& sz,
    true_type) // POD result type
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a2++, result++) {
      *result = ftor(a1, *a2, a3);
    }
  }

  // in_place_ftor(array, array)
  template <typename InPlaceBinaryFunctorType,
            typename ElementType1,
            typename ElementType2>
  inline void
  array_operation_in_place_binary(
    const InPlaceBinaryFunctorType& ftor,
    ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz)
  {
    ElementType1* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++, a2++) ftor(*a1, *a2);
  }

  // in_place_ftor(array, scalar)
  template <typename InPlaceBinaryFunctorType,
            typename ElementType>
  inline void
  array_operation_in_place_binary(
    const InPlaceBinaryFunctorType& ftor,
    ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz)
  {
    ElementType* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++) ftor(*a1, a2);
  }

  template <typename BooleanFunctorType,
            typename ElementType1,
            typename ElementType2>
  inline bool
  array_operation_reducing_boolean(
    const BooleanFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz)
  {
    const ElementType1* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++, a2++) {
      if (!ftor(*a1, *a2)) return false;
    }
    return true;
  }

  template <typename BooleanFunctorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean(
    const BooleanFunctorType& ftor,
    const ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz)
  {
    const ElementType* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++) {
      if (!ftor(*a1, a2)) return false;
    }
    return true;
  }

  template <typename BooleanFunctorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean(
    const BooleanFunctorType& ftor,
    const ElementType& a1,
    const ElementType* a2,
    const std::size_t& sz)
  {
    const ElementType* a2_end = a2 + sz;
    for(;a2 != a2_end; a2++) {
      if (!ftor(a1, *a2)) return false;
    }
    return true;
  }

  template <typename BooleanFunctorType,
            typename ElementType1,
            typename ElementType2>
  inline bool
  array_operation_reducing_boolean_not_equal_to(
    const BooleanFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz)
  {
    const ElementType1* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++, a2++) {
      if (ftor(*a1, *a2)) return true;
    }
    return false;
  }

  template <typename BooleanFunctorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_not_equal_to(
    const BooleanFunctorType& ftor,
    const ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz)
  {
    const ElementType* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++) {
      if (ftor(*a1, a2)) return true;
    }
    return false;
  }

  template <typename BooleanFunctorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_not_equal_to(
    const BooleanFunctorType& ftor,
    const ElementType& a1,
    const ElementType* a2,
    const std::size_t& sz)
  {
    const ElementType* a2_end = a2 + sz;
    for(;a2 != a2_end; a2++) {
      if (ftor(a1, *a2)) return true;
    }
    return false;
  }

  template <typename BooleanFunctorType,
            typename ElementType1,
            typename ElementType2>
  inline bool
  array_operation_reducing_boolean_greater_less(
    const BooleanFunctorType& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    const std::size_t& sz)
  {
    const ElementType1* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++, a2++) {
      if (ftor(*a1, *a2)) return true;
      if (ftor.reverse(*a2, *a1)) return false;
    }
    return false;
  }

  template <typename BooleanFunctorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_greater_less(
    const BooleanFunctorType& ftor,
    const ElementType* a1,
    const ElementType& a2,
    const std::size_t& sz)
  {
    const ElementType* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++) {
      if (ftor(*a1, a2)) return true;
      if (ftor.reverse(a2, *a1)) return false;
    }
    return false;
  }

  template <typename BooleanFunctorType,
            typename ElementType>
  inline bool
  array_operation_reducing_boolean_greater_less(
    const BooleanFunctorType& ftor,
    const ElementType& a1,
    const ElementType* a2,
    const std::size_t& sz)
  {
    const ElementType* a2_end = a2 + sz;
    for(;a2 != a2_end; a2++) {
      if (ftor(a1, *a2)) return true;
      if (ftor.reverse(*a2, a1)) return false;
    }
    return false;
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H
