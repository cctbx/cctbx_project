/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_C_INDEX_1D_CALCULATOR_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_C_INDEX_1D_CALCULATOR_H

#include <cstddef>

namespace scitbx { namespace af {

  template <std::size_t N>
  struct c_index_1d_calculator
  {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t
    get(ExtendArrayType const& e, IndexType const& i)
    {
      return c_index_1d_calculator<N-1>::get(e, i) * e[N-1] + i[N-1];
    }
  };

  template<>
  struct c_index_1d_calculator<1>
  {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t
    get(ExtendArrayType const& e, IndexType const& i)
    {
      return i[0];
    }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_C_INDEX_1D_CALCULATOR_H
