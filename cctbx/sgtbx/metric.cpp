// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/groups.h>

namespace cctbx { namespace sgtbx {

  bool SpaceGroup::isCompatibleMetricalMatrix(const af::double9& G,
                                              double tolerance) const
  {
    // for all R in the representative set of rotation matrices,
    // assert Transpose[R].G.R == G

    for (int iSMx = 1; iSMx < m_nSMx; iSMx++) {
      af::double9 R;
      int i;
      for(i=0;i<9;i++) R[i] = double(m_SMx[iSMx].Rpart()[i]);
      af::double9 RtGR = uctbx::getRtGR(G, R);
      for(i=0;i<9;i++) {
        double delta = RtGR[i] - G[i];
        if (delta < 0.) delta *= -1.;
        if (delta > tolerance) return false;
      }
    }
    return true;
  }

  void SpaceGroup::CheckMetricalMatrix(const af::double9& G,
                                       double tolerance) const
  {
    if (!isCompatibleMetricalMatrix(G, tolerance)) {
      throw error(
        "Unit cell is incompatible with symmetry operations.");
    }
  }

}} // namespace cctbx::sgtbx
