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

#include <stdexcept>

#include <cctbx/array_family/type_traits.h>

// XXX
#include <iostream>
#define CheckPoint std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl << std::flush

namespace cctbx { namespace af {

  inline
  void throw_range_error() {
    throw std::range_error("array_family");
  }

  struct reserve_flag {};
  struct no_initialization_flag {};

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

  // XXX use std::copy if compiler permits
  template <typename InputElementType,
            typename OutputElementType>
  OutputElementType*
  copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
    OutputElementType* p = result;
    while (first != last) *p++ = OutputElementType(*first++);
    return result;
  }

  // XXX use std::uninitialized_copy if compiler permits
  template <typename InputElementType,
            typename OutputElementType>
  OutputElementType*
  uninitialized_copy_typeconv(
    const InputElementType* first,
    const InputElementType* last,
    OutputElementType* result)
  {
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
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_MISC_H
