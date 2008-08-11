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
    get(ExtendArrayType const& /*e*/, IndexType const& i)
    {
      return i[0];
    }
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_C_INDEX_1D_CALCULATOR_H
