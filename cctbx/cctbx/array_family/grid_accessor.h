// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Feb 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_GRID_ACCESSOR_H
#define CCTBX_ARRAY_FAMILY_GRID_ACCESSOR_H

#include <cstddef>
#include <algorithm>
#include <cctbx/array_family/misc.h>
#include <cctbx/array_family/tiny_helpers.h>

// forward declaration
namespace cctbx { namespace af {
  template <typename ElementType, std::size_t N>
  class tiny;
}}

namespace cctbx { namespace af {

  template <std::size_t N>
  struct c_index_1d_calculator {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(const ExtendArrayType& e, const IndexType& i) {
      return c_index_1d_calculator<N-1>::get(e, i) * e[N-1] + i[N-1];
    }
  };

  template<>
  struct c_index_1d_calculator<1> {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(const ExtendArrayType& e, const IndexType& i) {
      return i[0];
    }
  };

  template <std::size_t N>
  struct fortran_index_1d_calculator {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(const ExtendArrayType& e, const IndexType& i) {
      return fortran_index_1d_calculator<N-1>::get(e, i)
        * e[e.size()-N] + i[e.size()-N];
    }
  };

  template<>
  struct fortran_index_1d_calculator<1> {
    template <typename ExtendArrayType, typename IndexType>
    static std::size_t get(const ExtendArrayType& e, const IndexType& i) {
      return i[e.size()-1];
    }
  };

  template <std::size_t N>
  struct compile_time_product {
    template <typename IndexType>
    static std::size_t get(const IndexType& i) {
      return i[N-1] * compile_time_product<N-1>::get(i);
    }
  };

  template <>
  struct compile_time_product<1> {
    template <typename IndexType>
    static std::size_t get(const IndexType& i) {
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
      typedef typename IndexType::value_type value_type;

      grid() { std::fill(this->begin(), this->end(), 0); }

      grid(const IndexType& n) : IndexType(n) {}

      CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(grid)

      static std::size_t nd() { return Nd; }

      void init_default() {
        for(std::size_t i=0;i<nd();i++) this->elems[i] = 0;
      }

      std::size_t size1d() const {
        return compile_time_product<Nd>::get(IndexType(*this));
      }

      std::size_t operator()(const IndexType& i) const {
        return Index1dCalculatorType::get(*this, i);
      }

      bool is_valid_index(const IndexType& i) const {
        for(std::size_t j=0;j<nd();j++) {
          if (i[j] < 0 || i[j] >= this->elems[j]) return false;
        }
        return true;
      }
  };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_GRID_ACCESSOR_H
