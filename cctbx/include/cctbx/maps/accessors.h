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
#include <cctbx/array_family/flex_grid_accessor.h>

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

      template <typename OtherArrayType>
      grid_p1(af::array_adaptor<OtherArrayType> const& a_a) : IndexType(a_a) {}

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

      IndexType i_1d_as_i_nd(const value_type& i_1d) const {
        return af::index_1d_as_c_index_nd(*this, i_1d);
      }
  };

  template <std::size_t Nd,
            typename Index1dCalculatorType = af::c_index_1d_calculator<Nd>,
            typename IndexType = af::tiny<long, Nd> >
  class padded_grid_p1
  {
    public:
      typedef IndexType index_type;
      typedef typename IndexType::value_type value_type;

      padded_grid_p1()
      {
        std::fill(n_logical_.begin(), n_logical_.end(), value_type(0));
        std::fill(n_physical_.begin(), n_physical_.end(), value_type(0));
      }

      padded_grid_p1(const IndexType& n_logical,
                     const IndexType& n_physical)
        : n_logical_(n_logical),
          n_physical_(n_physical)
      {}

      CCTBX_ARRAY_FAMILY_TINY_CONVENIENCE_CONSTRUCTORS(padded_grid_p1)

      static std::size_t nd() { return Nd; }

      std::size_t size1d() const {
        return af::compile_time_product<Nd>::get(n_physical_);
      }

      IndexType
      p1_i(const IndexType& i) const {
        IndexType result;
        for(std::size_t j=0;j<Nd;j++) {
          result[j] = i[j] % n_logical_[j];
          if (result[j] < 0) result[j] += n_logical_[j];
        }
        return result;
      }

      std::size_t operator()(const IndexType& i) const {
        return Index1dCalculatorType::get(n_physical_, p1_i(i));
      }

      const IndexType& n_logical() const { return n_logical_; }
      const IndexType& n_physical() const { return n_physical_; }

      af::flex_grid<> as_flex_grid() const // XXX move out
      {
        af::flex_grid_default_index_type grid;
        af::flex_grid_default_index_type layout;
        for(std::size_t i=0;i<n_logical_.size();i++) {
          grid.push_back(n_physical_[i]);
          layout.push_back(n_logical_[i]);
        }
        return af::flex_grid<>(grid).set_layout(layout);
      }

    protected:
      IndexType n_logical_;
      IndexType n_physical_;
  };

  //! 1-dimensional array index corresponding to Miller index element.
  template <typename IntegerType>
  inline IntegerType h_as_ih(IntegerType h, std::size_t N, bool PositiveOnly)
  {
    if (PositiveOnly) {
      if (0 > h || h >= N) return -1;
    }
    else {
      IntegerType m = (N - 1) / 2;
      if (-m > h || h > m) return -1;
      else if (h < 0) return h + N;
    }
    return h;
  }

  //! Miller index element corresponding to 1-dimensional array index.
  template <typename IntegerType>
  inline
  IntegerType
  ih_as_h(IntegerType ih, std::size_t n)
  {
    if (ih <= n/2) return ih;
    return ih - n;
  }

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_ACCESSORS_H
