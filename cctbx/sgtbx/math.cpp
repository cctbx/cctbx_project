// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 May 08 added: SmithNormalForm (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/fixes/cstdlib>
#include <cctbx/math/utils.h>
#include <cctbx/array_family/ref.h>
#include <cctbx/array_family/tiny.h>
#include <cctbx/sgtbx/math.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/basic/define_range.h>

namespace cctbx { namespace sgtbx {

  int iRowEchelonFormT(int *M, int mr, int mc, int *T, int tc)
  {
    // C++ version of RowEchelonFormT from the CrystGAP package
    //     (GAP Version 3.4.4).
    // B. Eick, F. Ga"hler and W. Nickel
    //     Computing Maximal Subgroups and Wyckoff Positions of Space Groups
    //     Acta Cryst. (1997). A53, 467 - 474

    af::ref<int, af::grid<2> > M2D(M, mr, mc);
    af::ref<int, af::grid<2> > T2D(T, mr, tc);
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
          int a = math::abs(M2D(k, j));
          if (a != 0 && a < math::abs(M2D(i, j))) {
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
    af::const_ref<int, af::grid<2> > M2D(M, nr, nc);
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

  void SolveHomRE1(const int REMx[3], const int IxIndep[2], af::int3 Sol[4])
  {
    // REMx must be in row echelon form with Rank 1.

    const int  TrialV[4][2] =
      {{ 1,  0 },
       { 0,  1 },
       { 1,  1 },
       { 1, -1 },
      };

    for (int iPV = 0; iPV < 4; iPV++) {
      int i;
      for(i=0;i<3;i++) Sol[iPV][i] = 0;
      for(i=0;i<2;i++) Sol[iPV][IxIndep[i]] = TrialV[iPV][i];
      cctbx_assert(iREBacksubst(REMx, 0, 2, 3, Sol[iPV].begin(), 0) > 0);
    }
  }

  int SmithNormalForm(int *M, int mr, int mc, int *P, int *Q)
  {
    int rr = mr;
    int rc = mc;

    if (P) MatrixLite::identity(P, mr);
    if (Q) MatrixLite::identity(Q, mc);

    for (;;)
    {
      rr = iRowEchelonFormT(M, rr, rc, P, mr);
      if (MatrixLite::isDiagonal(M, rr, rc)) break;
      MatrixLite::transpose(M, rr, rc);

      rc = iRowEchelonFormT(M, rc, rr, Q, mc);
      if (MatrixLite::isDiagonal(M, rc, rr)) break;
      MatrixLite::transpose(M, rc, rr);
    }

    if (Q) MatrixLite::transpose(Q, mc, mc);
    return rr;
  }

}} // namespace cctbx::sgtbx
