// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_TINY_H
#define CCTBX_ARRAY_FAMILY_TINY_H

#include <cctbx/array_family/tiny_base.h>
#include <cctbx/array_family/operator_traits.h>

namespace cctbx { namespace af {

  // Automatic allocation, fixed size, standard operators.
  template <typename ElementType, std::size_t N>
  class tiny : public tiny_base<ElementType, N>
  {
    public:
      CCTBX_ARRAY_FAMILY_TYPEDEFS

      tiny() {}

      CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(tiny)
      CCTBX_ARRAY_FAMILY_TINY_COPY_AND_ASSIGNMENT(tiny)
      CCTBX_ARRAY_FAMILY_TAKE_REF(this->elems, N)
  };

  template <typename ElementType, std::size_t N>
  inline
  tiny<
    typename unary_operator_traits<
      ElementType>::arithmetic, N>
  operator-(const tiny<ElementType, N>& a) {
    tiny<
      typename unary_operator_traits<
        ElementType>::arithmetic, N> result;
    for(std::size_t i=0;i<N;i++) result[i] = -a[i];
    return result;
  }

  template <typename ElementType, std::size_t N>
  inline
  tiny<bool, N>
  operator!(const tiny<ElementType, N>& a) {
    tiny<
      typename unary_operator_traits<
        ElementType>::logical, N> result;
    for(std::size_t i=0;i<N;i++) result[i] = !a[i];
    return result;
  }

  template <typename ElementTypeLHS,
            typename ElementTypeRHS,
            std::size_t N>
  inline
  tiny<
    typename binary_operator_traits<
      ElementTypeLHS, ElementTypeRHS>::arithmetic, N>
  operator+(
    const tiny<ElementTypeLHS, N>& lhs,
    const tiny<ElementTypeRHS, N>& rhs) {
    tiny<
      typename binary_operator_traits<
        ElementTypeLHS, ElementTypeRHS>::arithmetic, N> result;
    for(std::size_t i=0;i<N;i++) result[i] = lhs[i] + rhs[i];
    return result;
  }

  template <typename ElementTypeLHS,
            typename ElementTypeRHS,
            std::size_t N>
  inline
  tiny<
    typename binary_operator_traits<
      ElementTypeLHS, ElementTypeRHS>::logical, N>
  operator&&(
    const tiny<ElementTypeLHS, N>& lhs,
    const tiny<ElementTypeRHS, N>& rhs) {
    tiny<
      typename binary_operator_traits<
        ElementTypeLHS, ElementTypeRHS>::logical, N> result;
    for(std::size_t i=0;i<N;i++) result[i] = lhs[i] && rhs[i];
    return result;
  }

  template <typename ElementTypeLHS,
            typename ElementTypeRHS,
            std::size_t N>
  inline
  tiny<
    typename binary_operator_traits<
      ElementTypeLHS, ElementTypeRHS>::boolean, N>
  equal_to(
    const tiny<ElementTypeLHS, N>& lhs,
    const tiny<ElementTypeRHS, N>& rhs) {
    tiny<
      typename binary_operator_traits<
        ElementTypeLHS, ElementTypeRHS>::boolean, N> result;
    for(std::size_t i=0;i<N;i++) result[i] = lhs[i] == rhs[i];
    return result;
  }

  template <typename ElementTypeLHS,
            typename ElementTypeRHS,
            std::size_t N>
  inline
  typename binary_operator_traits<
    ElementTypeLHS, ElementTypeRHS>::boolean
  operator==(
    const tiny<ElementTypeLHS, N>& lhs,
    const tiny<ElementTypeRHS, N>& rhs) {
    for(std::size_t i=0;i<N;i++) if (!(lhs[i] == rhs[i])) {
      return ElementTypeLHS() != ElementTypeLHS();
    }
    return ElementTypeLHS() == ElementTypeLHS();
  }

}} // namespace cctbx::array_family

#endif // CCTBX_ARRAY_FAMILY_TINY_H
