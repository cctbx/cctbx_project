/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Feb: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_GRID_ACCESSOR_H
#define SCITBX_ARRAY_FAMILY_GRID_ACCESSOR_H

#include <cstddef>
#include <algorithm>
#include <scitbx/array_family/detail/tiny_helpers.h>
#include <scitbx/array_family/array_adaptor.h>

// forward declaration
namespace scitbx { namespace af {
  template <typename ElementType, std::size_t N>
  class tiny;
}}

namespace scitbx { namespace af {

  template <std::size_t N>
  struct c_index_1d_calculator {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(ExtendArrayType const& e, IndexType const& i) {
      return c_index_1d_calculator<N-1>::get(e, i) * e[N-1] + i[N-1];
    }
  };

  template<>
  struct c_index_1d_calculator<1> {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(ExtendArrayType const& e, IndexType const& i) {
      return i[0];
    }
  };

  template <std::size_t N>
  struct fortran_index_1d_calculator {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(ExtendArrayType const& e, IndexType const& i) {
      return fortran_index_1d_calculator<N-1>::get(e, i)
        * e[e.size()-N] + i[e.size()-N];
    }
  };

  template<>
  struct fortran_index_1d_calculator<1> {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(ExtendArrayType const& e, IndexType const& i) {
      return i[e.size()-1];
    }
  };

  template <std::size_t N>
  struct compile_time_product {
    template <typename IndexType>
    static std::size_t get(IndexType const& i) {
      return i[N-1] * compile_time_product<N-1>::get(i);
    }
  };

  template <>
  struct compile_time_product<1> {
    template <typename IndexType>
    static std::size_t get(IndexType const& i) {
      return i[0];
    }
  };

  template <std::size_t Nd,
            typename Index1dCalculatorType = c_index_1d_calculator<Nd>,
            typename IndexType = tiny<long, Nd> >
  class grid : public IndexType
  {
    public:
      typedef IndexType index_type;
      typedef typename IndexType::value_type index_value_type;
      typedef index_value_type value_type;

      grid() { std::fill(this->begin(), this->end(), value_type(0)); }

      grid(IndexType const& n) : IndexType(n) {}

      template <typename OtherArrayType>
      grid(af::array_adaptor<OtherArrayType> const& a_a) : IndexType(a_a) {}

      SCITBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(grid)

      static std::size_t nd() { return Nd; }

      std::size_t size_1d() const
      {
        return compile_time_product<Nd>::get(*this);
      }

      std::size_t operator()(IndexType const& i) const {
        return Index1dCalculatorType::get(*this, i);
      }

      bool is_valid_index(IndexType const& i) const {
        for(std::size_t j=0;j<nd();j++) {
          if (i[j] < 0 || i[j] >= this->elems[j]) return false;
        }
        return true;
      }
  };

  template <typename IndexNdType,
            typename Index1dType>
  inline
  IndexNdType
  index_1d_as_c_index_nd(IndexNdType const& n, Index1dType i_1d)
  {
    IndexNdType i_nd;
    for(std::size_t j=n.size()-1;j;j--) {
      i_nd[j] = i_1d % n[j];
      i_1d /= n[j];
    }
    i_nd[0] = i_1d;
    return i_nd;
  }

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_GRID_ACCESSOR_H
