// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: moved parts from ref.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_MISC_H
#define CCTBX_ARRAY_FAMILY_MISC_H

#include <cctbx/array_family/type_traits.h>

namespace cctbx { namespace af {

  struct reserve_flag {};

  namespace detail {

    template <class ElementType>
    inline
    void destroy_array_element(ElementType* elem, false_type) {
      elem->~ElementType();
    }

    template <class ElementType>
    inline
    void destroy_array_element(ElementType* elem, true_type) {
    }

    template <class ElementType>
    inline
    void destroy_array_elements(ElementType* first, ElementType* last,
                                false_type) {
      while (first != last) {
        first->~ElementType();
        ++first;
      }
    }

    template <class ElementType>
    inline
    void destroy_array_elements(ElementType* first, ElementType* last,
                                true_type) {
    }

  } // namespace detail

#if !defined(__GNUC__) \
    || ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 0)))
#define CCTBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
#endif

  template <typename InputElementType,
            typename OutputElementType>
  inline
  OutputElementType*
  copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
#ifdef CCTBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
    return std::copy(first, last, result);
#else
    OutputElementType* p = result;
    while (first != last) *p++ = OutputElementType(*first++);
    return result;
#endif
  }

  template <typename InputElementType,
            typename OutputElementType>
#ifdef CCTBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
  inline
#endif
  OutputElementType*
  uninitialized_copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
#ifdef CCTBX_ARRAY_FAMILY_STD_COPY_PERMITS_TYPECONV
    return std::uninitialized_copy(first, last, result);
#else
    OutputElementType* p = result;
    try {
      for (; first != last; p++, first++) {
        new (p) OutputElementType(*first);
      }
    }
    catch (...) {
      detail::destroy_array_elements(result, p, false_type());
      throw;
    }
    return result;
#endif
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_MISC_H
