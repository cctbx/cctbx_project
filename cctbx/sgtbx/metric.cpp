// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/groups.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  bool SgOps::isCompatibleMetricalMatrix(const uctbx::Mx33& G,
                                         double tolerance) const
  {
    // for all R in the representative set of rotation matrices,
    // assert Transpose[R].G.R == G

    for (int iSMx = 1; iSMx < m_nSMx; iSMx++) {
      uctbx::Mx33 R;
      int i;
      for(i=0;i<9;i++) R[i] = static_cast<double>(m_SMx[iSMx].Rpart()[i]);
      uctbx::Mx33 RtGR = uctbx::getRtGR(G, R);
      for(i=0;i<9;i++) {
        double delta = RtGR[i] - G[i];
        if (delta < 0.) delta *= -1.;
        if (delta > tolerance) return false;
      }
    }
    return true;
  }

  void SgOps::CheckMetricalMatrix(const uctbx::Mx33& G,
                                  double tolerance) const
  {
    if (! isCompatibleMetricalMatrix(G, tolerance)) {
      throw error(
        "Unit cell is incompatible with symmetry operations.");
    }
  }

} // namespace sgtbx
