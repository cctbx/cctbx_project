// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPS_SYM_FLAGS_H
#define CCTBX_MAPS_SYM_FLAGS_H

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
      result[i] = result[i] * M.TBF() + M.Tpart()[i] * M.RBF() * N[i];
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

  template <typename VecRefNdType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(const VecRefNdType& p1_flags,
             const sgtbx::SpaceGroup& SgOps,
             const IndexTupleType& pivot)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_flags.dim()(pivot);
    for(int iSMx=1;iSMx<SgOps.OrderZ();iSMx++) {
      sgtbx::RTMx M = SgOps(iSMx);
      tagged_value<IndexTupleType>
      sym_equiv_point = multiply(p1_flags.dim(), M, pivot);
      if (sym_equiv_point.tag) {
        std::size_t i1d_sep = p1_flags.dim()(sym_equiv_point.value);
        if (i1d_sep != i1d_pivot) p1_flags[i1d_sep] = i1d_pivot;
      }
      else {
        grid_misses++;
      }
    }
    return grid_misses;
  }

  template <typename VecRefNdType,
            typename GridSsType,
            typename IndexTupleType>
  std::size_t
  mark_orbit(const VecRefNdType& p1_flags,
             const GridSsType& grid_ss,
             const IndexTupleType& pivot)
  {
    std::size_t grid_misses = 0;
    std::size_t i1d_pivot = p1_flags.dim()(pivot);
    fixcap_vector<int, 3> moduli(grid_ss.size());
    for(int i=0;i<grid_ss.size();i++) moduli[i] = grid_ss[i].M;
    nested_loop<fixcap_vector<int, 3> > loop(moduli);
    for (fixcap_vector<int, 3> f = loop(); !loop.over(); f = loop.next()) {
      tagged_value<IndexTupleType>
      sym_equiv_point = add(p1_flags.dim(), grid_ss, pivot, f);
      if (sym_equiv_point.tag) {
        std::size_t i1d_sep = p1_flags.dim()(sym_equiv_point.value);
        if (i1d_sep != i1d_pivot) p1_flags[i1d_sep] = i1d_pivot;
      }
      else {
        grid_misses++;
      }
    }
    return grid_misses;
  }

  template <typename VecRefNdType,
            typename SymmetryType>
  std::size_t
  mark_orbits(const VecRefNdType& p1_flags,
              const SymmetryType& symmetry)
  {
    std::size_t grid_misses = 0;
    nested_loop<int3> loop(p1_flags.dim());
    for (int3 pivot = loop(); !loop.over(); pivot = loop.next()) {
      if (p1_flags(pivot) == -1) {
        grid_misses += mark_orbit(p1_flags, symmetry, pivot);
      }
    }
    return grid_misses;
  }

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_SYM_FLAGS_H
