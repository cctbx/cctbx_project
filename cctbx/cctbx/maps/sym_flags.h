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

  template <typename DimTupleType,
            typename IndexTupleType>
  boost::array<int, 3>
  multiply(const DimTupleType& N,
           const sgtbx::RTMx& M,
           const IndexTupleType& X)
  {
    boost::array<int, 3> result = M.Rpart() * X;
    for(int i=0;i<3;i++) {
      result[i] = result[i] * M.TBF() + M.Tpart()[i] * M.RBF() * N[i];
      if (result[i] % (M.RBF() * M.TBF())) {
        throw error("Grid is not compatible with symmetry.");
      }
      result[i] /= (M.RBF() * M.TBF());
    }
    return result;
  }

  template <typename VecRefNdType,
            typename IndexTupleType>
  void
  mark_orbit(const sgtbx::SpaceGroup& SgOps,
             const VecRefNdType& p1_flags,
             const IndexTupleType& pivot)
  {
    std::size_t i1d_pivot = p1_flags.dim()(pivot);
    for(int iSMx=1;iSMx<SgOps.OrderZ();iSMx++) {
      sgtbx::RTMx M = SgOps(iSMx);
      boost::array<int, 3>
      sym_equiv_point = multiply(p1_flags.dim(), M, pivot);
      std::size_t i1d_sep = p1_flags.dim()(sym_equiv_point);
      if (i1d_sep != i1d_pivot) p1_flags[i1d_sep] = i1d_pivot;
    }
  }

}} // namespace cctbx::maps

#endif // CCTBX_MAPS_SYM_FLAGS_H
