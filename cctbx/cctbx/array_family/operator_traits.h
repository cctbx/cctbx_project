/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (Ralf W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_OPERATOR_TRAITS_H
#define CCTBX_ARRAY_FAMILY_OPERATOR_TRAITS_H

#include <cctbx/array_family/operator_traits_builtin.h>
#include <cctbx/array_family/flagged_value.h>

namespace cctbx { namespace af {

  template<typename ValueTypeLhs, typename ValueTypeRhs>
  struct binary_operator_traits<
    flagged_value<ValueTypeLhs>,
    flagged_value<ValueTypeRhs> > {
    typedef binary_operator_traits<ValueTypeLhs, ValueTypeRhs> val_traits;
    typedef flagged_value<typename val_traits::arithmetic> arithmetic;
    typedef flagged_value<typename val_traits::logical> logical;
    typedef flagged_value<typename val_traits::boolean> boolean;
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_OPERATOR_TRAITS_H
