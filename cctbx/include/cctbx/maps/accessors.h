// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_ACCESSORS_H
#define CCTBX_MAPS_ACCESSORS_H

#include <cctbx/array_family/grid_accessor.h>

namespace cctbx {

//! Algorithms for the handling of 3-dimensional %maps.
namespace maps {


  template <std::size_t Nd,
            typename Index1dCalculatorType = af::c_index_1d_calculator<Nd>,
            typename IndexType = af::tiny<long, Nd> >
  class grid_p1 : public IndexType
  {
    public:
      typedef IndexType index_type;
      typedef typename IndexType::value_type value_type;

      grid_p1() { std::fill(this->begin(), this->end(), value_type(0)); }

      grid_p1(const IndexType& n) : IndexType(n) {}

      CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(grid_p1)

      static std::size_t nd() { return Nd; }

      std::size_t size1d() const {
        return af::compile_time_product<Nd>::get(IndexType(*this));
      }

      IndexType
      p1_i(const IndexType& i) const {
        IndexType result;
        for(std::size_t j=0;j<Nd;j++) {
          result[j] = i[j] % this->elems[j];
          if (result[j] < 0) result[j] += this->elems[j];
        }
        return result;
      }

      std::size_t operator()(const IndexType& i) const {
        return Index1dCalculatorType::get(*this, p1_i(i));
      }
  };

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_ACCESSORS_H
