#ifndef SCITBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H
#define SCITBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H

#include <scitbx/array_family/detail/misc.h>

namespace scitbx { namespace af {

  // functor(array), non-POD result type
  template <typename FunctorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_a(
    FunctorType const& ftor,
    const ElementType* a,
    ElementTypeResult* result,
    std::size_t const& sz,
    false_type)
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a++, result++) {
        new (result) ElementTypeResult(ftor(*a));
      }
    }
    catch (...) {
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array), POD result type
  template <typename FunctorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_a(
    FunctorType const& ftor,
    const ElementType* a,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a++, result++) {
      *result = ftor(*a);
    }
  }

  // functor(array, array), non-POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_a_a(
    FunctorType const& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    std::size_t const& sz,
    false_type)
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, a2++, result++) {
        new (result) ElementTypeResult(ftor(*a1, *a2));
      }
    }
    catch (...) {
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array, array), POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_a_a(
    FunctorType const& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, a2++, result++) {
      *result = ftor(*a1, *a2);
    }
  }

  // functor(array, scalar), non-POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_a_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    ElementType2 const& a2,
    ElementTypeResult* result,
    std::size_t const& sz,
    false_type)
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, result++) {
        new (result) ElementTypeResult(ftor(*a1, a2));
      }
    }
    catch (...) {
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array, scalar), POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_a_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    ElementType2 const& a2,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, result++) {
      *result = ftor(*a1, a2);
    }
  }

  // functor(scalar, array), non-POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  void // not inline
  array_operation_s_a(
    FunctorType const& ftor,
    ElementType1 const& a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    std::size_t const& sz,
    false_type)
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a2++, result++) {
        new (result) ElementTypeResult(ftor(a1, *a2));
      }
    }
    catch (...) {
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(scalar, array), POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementTypeResult>
  inline void
  array_operation_s_a(
    FunctorType const& ftor,
    ElementType1 const& a1,
    const ElementType2* a2,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a2++, result++) {
      *result = ftor(a1, *a2);
    }
  }

  // functor(array, array, scalar), non-POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  void // not inline
  array_operation_a_a_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementType3 const& a3,
    ElementTypeResult* result,
    std::size_t const& sz,
    false_type)
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, a2++, result++) {
        new (result) ElementTypeResult(ftor(*a1, *a2, a3));
      }
    }
    catch (...) {
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array, array, scalar), POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  inline void
  array_operation_a_a_s(
    FunctorType const& ftor,
    const ElementType1* a1,
    const ElementType2* a2,
    ElementType3 const& a3,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, a2++, result++) {
      *result = ftor(*a1, *a2, a3);
    }
  }

  // functor(array, scalar, scalar), non-POD result type
  template <typename FunctorType,
            typename ElementType,
            typename ElementTypeResult>
  void // not inline
  array_operation_a_s_s(
    FunctorType const& ftor,
    const ElementType* a1,
    ElementType const& a2,
    ElementType const& a3,
    ElementTypeResult* result,
    std::size_t const& sz,
    false_type)
  {
    ElementTypeResult* result_start = result;
    try {
      ElementTypeResult* result_end = result + sz;
      for(;result != result_end; a1++, result++) {
        new (result) ElementTypeResult(ftor(*a1, a2, a3));
      }
    }
    catch (...) {
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(array, scalar, scalar), POD result type
  template <typename FunctorType,
            typename ElementType,
            typename ElementTypeResult>
  inline void
  array_operation_a_s_s(
    FunctorType const& ftor,
    const ElementType* a1,
    ElementType const& a2,
    ElementType const& a3,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a1++, result++) {
      *result = ftor(*a1, a2, a3);
    }
  }

  // functor(scalar, array, scalar), non-POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  void // not inline
  array_operation_s_a_s(
    FunctorType const& ftor,
    ElementType1 const& a1,
    const ElementType2* a2,
    ElementType3 const& a3,
    ElementTypeResult* result,
    std::size_t const& sz,
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
      scitbx::af::detail::destroy_array_elements(result_start, result,
        false_type());
      throw;
    }
  }

  // functor(scalar, array, scalar), POD result type
  template <typename FunctorType,
            typename ElementType1,
            typename ElementType2,
            typename ElementType3,
            typename ElementTypeResult>
  inline void
  array_operation_s_a_s(
    FunctorType const& ftor,
    ElementType1 const& a1,
    const ElementType2* a2,
    ElementType3 const& a3,
    ElementTypeResult* result,
    std::size_t const& sz,
    true_type)
  {
    ElementTypeResult* result_end = result + sz;
    for(;result != result_end; a2++, result++) {
      *result = ftor(a1, *a2, a3);
    }
  }

  // functor in_place(array, array)
  template <typename InPlaceFunctorType,
            typename ElementType1,
            typename ElementType2>
  inline void
  array_operation_in_place_a_a(
    InPlaceFunctorType const& ftor,
    ElementType1* a1,
    const ElementType2* a2,
    std::size_t const& sz)
  {
    ElementType1* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++, a2++) ftor(*a1, *a2);
  }

  // functor in_place(array, scalar)
  template <typename InPlaceFunctorType,
            typename ElementType>
  inline void
  array_operation_in_place_a_s(
    InPlaceFunctorType const& ftor,
    ElementType* a1,
    ElementType const& a2,
    std::size_t const& sz)
  {
    ElementType* a1_end = a1 + sz;
    for(;a1 != a1_end; a1++) ftor(*a1, a2);
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_GENERIC_ARRAY_OPERATORS_H
