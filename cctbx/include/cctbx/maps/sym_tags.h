// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_SYM_TAGS_H
#define CCTBX_MAPS_SYM_TAGS_H

#include <cctbx/array_family/tiny_types.h>
#include <cctbx/array_family/small.h>
#include <cctbx/array_family/loops.h>
#include <cctbx/math/linear_regression.h>
#include <cctbx/sgtbx/groups.h>

namespace cctbx { namespace maps {

  template <typename ValueType, typename TagType = bool>
  struct tagged_value
  {
    tagged_value() {}
    tagged_value(const ValueType& v)
      : value(v)
    {}
    tagged_value(const ValueType& v, const TagType& t)
      : value(v), tag(t)
    {}
    ValueType value;
    TagType tag;
  };

  template <typename DimTupleType,
            typename IndexTupleType>
  inline tagged_value<IndexTupleType>
  multiply(const DimTupleType& N,
           const sgtbx::RTMx& M,
           const IndexTupleType& X)
  {
    IndexTupleType result = M.Rpart() * X;
    for(int i=0;i<3;i++) {
      result[i] = result[i] * M.TBF() + M.Tpart().vec()[i] * M.RBF() * N[i];
      if (result[i] % (M.RBF() * M.TBF())) {
        return tagged_value<IndexTupleType>(result, false);
      }
      result[i] /= (M.RBF() * M.TBF());
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename DimensionType,
            typename GridSsType,
            typename IndexTupleType,
            typename FactorTupleType>
  inline tagged_value<IndexTupleType>
  add(const DimensionType& dim,
      const GridSsType& grid_ss,
      const IndexTupleType& pivot,
      const FactorTupleType& f)
  {
    IndexTupleType result = pivot;
    for(int i_ss=0;i_ss<grid_ss.size();i_ss++) {
      for(int i=0;i<3;i++) {
        int s = dim[i] * grid_ss[i_ss].V[i] * f[i_ss];
        if (s % grid_ss[i_ss].M) {
          return tagged_value<IndexTupleType>(result, false);
        }
        result[i] += s / grid_ss[i_ss].M;
      }
    }
    return tagged_value<IndexTupleType>(result, true);
  }

  template <typename TagVersaType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(TagVersaType& p1_tags,
             const sgtbx::SpaceGroup& SgOps,
             const IndexTupleType& pivot)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    for(int iSMx=1;iSMx<SgOps.OrderZ();iSMx++) {
      sgtbx::RTMx M = SgOps(iSMx);
      tagged_value<IndexTupleType>
      sym_equiv_point = multiply(p1_tags.accessor(), M, pivot);
      if (sym_equiv_point.tag) {
        std::size_t i1d_sep = p1_tags.accessor()(sym_equiv_point.value);
        while (p1_tags[i1d_sep] != -1) i1d_sep = p1_tags[i1d_sep];
        if (i1d_sep != i1d_pivot) p1_tags[i1d_sep] = i1d_pivot;
      }
      else {
        grid_misses++;
      }
    }
    return grid_misses;
  }

  template <typename TagVersaType,
            typename GridSsType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(TagVersaType& p1_tags,
             const GridSsType& grid_ss,
             const IndexTupleType& pivot)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_tags.accessor()(pivot);
    af::small<int, 3> moduli(grid_ss.size());
    for(int i=0;i<grid_ss.size();i++) moduli[i] = grid_ss[i].M;
    af::nested_loop<af::small<int, 3> > loop(moduli);
    for (const af::small<int, 3>& f = loop(); !loop.over(); loop.incr()) {
      tagged_value<IndexTupleType>
      sym_equiv_point = add(p1_tags.accessor(), grid_ss, pivot, f);
      if (sym_equiv_point.tag) {
        std::size_t i1d_sep = p1_tags.accessor()(sym_equiv_point.value);
        while (p1_tags[i1d_sep] != -1) i1d_sep = p1_tags[i1d_sep];
        if (i1d_sep != i1d_pivot) p1_tags[i1d_sep] = i1d_pivot;
      }
      else {
        grid_misses++;
      }
    }
    return grid_misses;
  }

  template <typename TagVersaType,
            typename SymmetryType>
  std::size_t
  mark_orbits(TagVersaType& p1_tags,
              const SymmetryType& symmetry)
  {
    std::size_t grid_misses = 0;
    af::nested_loop<af::int3> loop(p1_tags.accessor());
    for (const af::int3& pivot = loop(); !loop.over(); loop.incr()) {
      if (p1_tags(pivot) == -1) {
        grid_misses += mark_orbit(p1_tags, symmetry, pivot);
      }
    }
    return grid_misses;
  }

  template <typename TagArrayType>
  std::size_t
  optimize_tags(TagArrayType tags)
  {
    std::size_t n_independent = 0;
    for(std::size_t i=0;i<tags.size();i++) {
      if (tags[i] < 0) {
        n_independent++;
      }
      else {
        std::size_t j = tags[i];
        while (tags[j] >= 0) j = tags[j];
        tags[i] = j;
      }
    }
    return n_independent;
  }

  template <class DataArrayType>
  struct verify_grid_tags
    : math::linear_regression_core<typename DataArrayType::size_type,
                                   typename DataArrayType::value_type>
  {
    typedef typename DataArrayType::size_type size_type;
    typedef typename DataArrayType::value_type value_type;

    verify_grid_tags() {}
    template <typename TagArrayType>
    verify_grid_tags(const DataArrayType& data,
                     const TagArrayType& tags,
                     const value_type& epsilon = 1.e-6)
    {
      cctbx_assert(data.size() <= tags.size());
      n_dependent = 0;
      size_type i;
      for(i=0;i<data.size();i++) {
        if (tags[i] >= 0) break;
      }
      if (i == data.size()) {
        reset();
        return;
      }
      n_dependent = 1;
      value_type x = data[i];
      value_type y = data[tags[i]];
      value_type min_x = x;
      value_type max_x = x;
      value_type min_y = y;
      value_type max_y = y;
      value_type sum_x = x;
      value_type sum_x2 = x * x;
      value_type sum_y = y;
      value_type sum_y2 = y * y;
      value_type sum_xy = x * y;
      for(i++;i<data.size();i++) {
        if (tags[i] < 0) continue;
        n_dependent++;
        x = data[i];
        y = data[tags[i]];
        if (min_x > x) min_x = x;
        if (max_x < x) max_x = x;
        if (min_y > y) min_y = y;
        if (max_y < y) max_y = y;
        sum_x += x;
        sum_x2 += x * x;
        sum_y += y;
        sum_y2 += y * y;
        sum_xy += x * y;
      }
      set(n_dependent,
          min_x, max_x, min_y, max_y, sum_x, sum_x2, sum_y, sum_y2, sum_xy,
          epsilon);
    }
    size_type n_dependent;
  };

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_SYM_TAGS_H
