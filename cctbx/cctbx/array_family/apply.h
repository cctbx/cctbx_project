// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_APPLY_H
#define CCTBX_ARRAY_FAMILY_APPLY_H

#include <cctbx/basic/meta.h>

namespace cctbx { namespace af {

  template <typename ArrayType, typename OtherElementType>
  struct change_array_element_type {
    typedef void array_type;
  };

  template <typename UnaryOperation,
            typename ArgumentArrayType,
            typename ResultElementType>
  typename
  change_array_element_type<ArgumentArrayType, ResultElementType>::array_type
  apply(UnaryOperation op,
        const ArgumentArrayType& arguments,
        type_holder<ResultElementType> result_type_holder)
  {
    typename
    change_array_element_type<ArgumentArrayType, ResultElementType>::array_type
    result(arguments.size());
    for (std::size_t i=0;i<arguments.size();i++) {
      result[i] = op(arguments[i]);
    }
    return result;
  }

  template <typename UnaryOperation,
            typename ArgumentArrayType>
  typename
  change_array_element_type<
    ArgumentArrayType,
    typename UnaryOperation::result_type>::array_type
  apply(UnaryOperation op,
        const ArgumentArrayType& arguments)
  {
    return apply(op, arguments,
                 type_holder<typename UnaryOperation::result_type>());
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_APPLY_H
