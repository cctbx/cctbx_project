// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/math.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  int iRowEchelonFormT(int *M, int mr, int mc, int *T, int tc)
  {
    // C++ version of RowEchelonFormT from the CrystGAP package
    //     (GAP Version 3.4.4).
    // B. Eick, F. Ga"hler and W. Nickel
    //     Computing Maximal Subgroups and Wyckoff Positions of Space Groups
    //     Acta Cryst. (1997). A53, 467 - 474

    Array2DAdaptor<int> M2D(M, mc);
    Array2DAdaptor<int> T2D(T, tc);
    int i, j;
    for (i = j = 0; i < mr && j < mc;) {
      int k = i; while (k < mr && M2D(k, j) == 0) k++;
      if (k == mr)
        j++;
      else {
        if (i != k) {
                 swap(&M2D(i, 0), &M2D(k, 0), mc);
          if (T) swap(&T2D(i, 0), &T2D(k, 0), tc);
        }
        for (k++; k < mr; k++) {
          int a = abs(M2D(k, j));
          if (a != 0 && a < abs(M2D(i, j))) {
                   swap(&M2D(i, 0), &M2D(k, 0), mc);
            if (T) swap(&T2D(i, 0), &T2D(k, 0), tc);
          }
        }
        if (M2D(i, j) < 0) {
          int ic;
                 for(ic=0;ic<mc;ic++) M2D(i, ic) *= -1;
          if (T) for(ic=0;ic<tc;ic++) T2D(i, ic) *= -1;
        }
        int Cleared = 1;
        for (k = i + 1; k < mr; k++) {
          int a = M2D(k, j) / M2D(i, j);
          if (a != 0) {
            int ic;
                   for(ic=0;ic<mc;ic++) M2D(k, ic) -= a * M2D(i, ic);
            if (T) for(ic=0;ic<tc;ic++) T2D(k, ic) -= a * T2D(i, ic);
          }
          if (M2D(k, j) != 0) Cleared = 0;
        }
        if (Cleared) { i++; j++; }
      }
    }
    return i;
  }

  int iREBacksubst(const int *M, const int *V,
                   const int nr, const int nc,
                   int *Sol, int *FlagIndep)
  {
    Array2DAdaptor<const int> M2D(M, nc);
    if (FlagIndep) {
      for (int ic = 0; ic < nc; ic++) FlagIndep[ic] = 1;
    }
    int d = 1;
    for (int ir = nr - 1; ir >= 0; ir--)
    {
      int ic;
      for (ic = 0; ic < nc; ic++) {
        if (M2D(ir, ic)) goto Set_Sol_ic;
      }
      if (V && V[ir] != 0) return 0;
      continue;

      Set_Sol_ic:
      if (FlagIndep) FlagIndep[ic] = 0;
      if (Sol) {
                  int icp = ic + 1;
        int nv = nc - icp;
        if (nv) {
          MatrixLite::multiply<int>(&M2D(ir, icp), &Sol[icp], 1, nv, 1,
                                    &Sol[ic]);
          Sol[ic] *= -1;
        }
        else {
          Sol[ic] = 0;
        }
        if (V) Sol[ic] += d * V[ir];
                      int m = M2D(ir, ic);
        int f = gcd(Sol[ic], m);
        if (m < 0) f *= -1;
        Sol[ic] /= f;
        f = m / f;
        if (f != 1) {
          int jc;
          for(jc=0;jc<nc;jc++) if (jc != ic) Sol[jc] *= f;
          d *= f;
        }
      }
    }
    return d;
  }

  int iRESetIxIndep(const int *REMx, int nr, int nc, int *IxIndep, int mIndep)
  {
    int FlagIndep[6];
    cctbx_assert(nc <= sizeof FlagIndep / sizeof (*FlagIndep));
    cctbx_assert(iREBacksubst(REMx, 0, nr, nc, 0, FlagIndep) >= 1);
    int nIndep = 0;
    int ic;
    for(ic=0;ic<nc;ic++) {
      if (FlagIndep[ic]) {
        if (nIndep == mIndep) return -1;
        IxIndep[nIndep++] = ic;
      }
    }
    return nIndep;
  }

} // namespace sgtbx
