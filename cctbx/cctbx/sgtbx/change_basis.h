// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_CHANGE_BASIS_H
#define CCTBX_SGTBX_CHANGE_BASIS_H

#include <cctbx/sgtbx/matrix.h>

namespace sgtbx {

  class ChOfBasisOp {
    private:
      RTMx Mx;
      RTMx InvMx;
    public:
      inline ChOfBasisOp(const RTMx& M, const RTMx& InvM)
        : Mx(M), InvMx(InvM) { }
      explicit inline ChOfBasisOp(const RTMx& M)
        : Mx(M), InvMx(M.inverse()) { }
      explicit inline ChOfBasisOp(int RBF = CRBF)
        : Mx(RBF), InvMx(RBF) { }
      explicit inline ChOfBasisOp(int RBF = CRBF, int TBF = CTBF)
        : Mx(RBF, TBF), InvMx(RBF, TBF) { }
      inline bool isValid() const {
        return Mx.isValid() && InvMx.isValid();
      }
      inline const RTMx& M() const { return Mx; }
      inline const RTMx& InvM() const { return InvMx; }
      inline ChOfBasisOp swap() const { return ChOfBasisOp(InvMx, Mx); }
      // Mx * M * InvM
      RTMx operator()(const RTMx& RT) const;
      // Mx * (I|T) * InvM
      TrVec operator()(const TrVec& T, int SignI) const;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_CHANGE_BASIS_H
