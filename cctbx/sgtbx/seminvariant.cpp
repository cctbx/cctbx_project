// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 02: start port of sglite/sgss.c (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/reference.h>

#if 1
#include <iostream>
using std::cout;
using std::endl;
#define CheckPoint cout << __FILE__ << "(" << __LINE__ << ")" << endl
#endif

namespace sgtbx {

#ifdef JUNK
  const int mDiscrGr = 8;

  struct T_DiscrGr
    {
      int  P[3]; // primitive cell
      int  Z[3]; // centred cell
    };

  int Is_ss(const T_ssVM *ssVM, int n_ssVM, int h, int k, int l)
  {
    int  i_ssVM, u;

    range1(i_ssVM, n_ssVM)
    {
      u =  ssVM[i_ssVM].V[0] * h;
      u += ssVM[i_ssVM].V[1] * k;
      u += ssVM[i_ssVM].V[2] * l;

      if (ssVM[i_ssVM].M) {
        if (u % ssVM[i_ssVM].M) return 0; }
      else {
        if (u)                  return 0; }
    }

    return 1;
  }


  void Set_uvw(const T_ssVM *ssVM, int n_ssVM, int h, int k, int l, int uvw[3])
  {
    int  i_ssVM, u;

    range1(i_ssVM, n_ssVM)
    {
      u =  ssVM[i_ssVM].V[0] * h;
      u += ssVM[i_ssVM].V[1] * k;
      u += ssVM[i_ssVM].V[2] * l;

      if (ssVM[i_ssVM].M) u %= ssVM[i_ssVM].M;

      uvw[i_ssVM] = u;
    }
  }
#endif

  namespace detail {

    struct AnyGenerators {
      AnyGenerators(const SgOps& sgo);
      void setPrimitive();
      ChOfBasisOp Z2POp;
      TrVec ZInvT;
      TrVec PInvT;
      int nGen;
      RTMx ZGen[2];
      RTMx PGen[2];
    };

    AnyGenerators::AnyGenerators(const SgOps& sgo)
      : nGen(0)
    {
      using namespace tables::CrystalSystem;

      Z2POp = sgo.getZ2POp();

      ZInvT = sgo.InvT(true);
      PInvT = TrVec(0);
      int i;
      for(i=0;i<2;i++) ZGen[i] = RTMx(0, 0);
      for(i=0;i<2;i++) PGen[i] = RTMx(0, 0);

      int PrincipalProperOrder = 0;

      tables::MatrixGroup::Code PG_MGC = sgo.getPointGroupType();
      switch (PG_MGC.CrystalSystem())
      {
        case Triclinic:
          break;

        case Monoclinic:
          ZGen[0] = sgo[1];
          nGen = 1;
          break;

        case Orthorhombic:
          ZGen[0] = sgo[1];
          ZGen[1] = sgo[2];
          nGen = 2;
          break;

        case Tetragonal:
                                     PrincipalProperOrder = 4;
        case Trigonal:
          if (!PrincipalProperOrder) PrincipalProperOrder = 3;
        case Hexagonal:
          if (!PrincipalProperOrder) PrincipalProperOrder = 6;
          {
            RotMxInfo PrincipalRI;
            for(i = 1; i < sgo.nSMx(); i++) {
              PrincipalRI = sgo[i].Rpart().getInfo();
              if (std::abs(PrincipalRI.Rtype()) == PrincipalProperOrder) {
                if (PrincipalRI.SenseOfRotation() > 0) {
                  ZGen[0] = sgo[i];
                  nGen++;
                  break;
                }
              }
            }
            cctbx_assert(nGen == 1);
            int iPrincipal = i;
            for(i = 1; i < sgo.nSMx(); i++) {
              if (i == iPrincipal) continue;
              RotMxInfo RI = sgo[i].Rpart().getInfo();
              if (std::abs(RI.Rtype()) == 2) {
                if (PrincipalRI.EV() != RI.EV()) {
                  ZGen[1] = sgo[i];
                  nGen++;
                  break;
                }
              }
            }
          }
          break;

        case Cubic:
          for(i = 1; i < sgo.nSMx(); i++) {
            RotMxInfo RI = sgo[i].Rpart().getInfo();
            if      (std::abs(RI.Rtype()) == 3) {
              if (RI.SenseOfRotation() > 0) {
                if (!ZGen[0].isValid()) {
                  ZGen[0] = sgo[i];
                  nGen++;
                  if (nGen == 2) break;
                }
              }
            }
            else if (std::abs(RI.Rtype()) == sgo.nSMx() / 6) {
              if (RI.SenseOfRotation() >= 0) {
                if (!ZGen[1].isValid()) {
                  ZGen[1] = sgo[i];
                  nGen++;
                  if (nGen == 2) break;
                }
              }
            }
          }
          cctbx_assert(nGen == 2);
          break;

        default:
          throw cctbx_internal_error();
      }
    }

    void AnyGenerators::setPrimitive()
    {
      for (int i = 0; i < nGen; i++) {
        PGen[i] = Z2POp(ZGen[i]).modPositive();
      }
      if (ZInvT.isValid()) {
        PInvT = Z2POp(ZInvT, -1).modPositive();
      }
    }

  } // namespace detail


#ifdef JUNK
  static int CmpOLen2(const int a[3], const int b[3])
  {
    int OLen2a, OLen2b, i;

    OLen2a = 0; rangei(3) OLen2a += a[i] * a[i];
    OLen2b = 0; rangei(3) OLen2b += b[i] * b[i];

    if (OLen2a < OLen2b) return -1;
    if (OLen2a > OLen2b) return  1;

    return CmpiVect(a, b, 3);
  }


  static int ConstructGenRmI(const T_SgOps *SgOps, const T_RTMx Z2PCBMx[2],
                             const int IxGen[2], int nGen,
                             int GenRmI[3 * 3 * 3])
  {
    int  nrGenRmI;
    int  iRmI, iGen, i;

        nrGenRmI = (SgOps->fInv - 1 + nGen) * 3;
    if (nrGenRmI > 9) return IE(-1);

    iRmI = 0;

    if (SgOps->fInv == 2) {
      SetRminusI(SgOps->SMx[0].s.R, GenRmI, 1);
      iRmI++;
    }

    if (Z2PCBMx == NULL) {
      range1(iGen, nGen) {
        SetRminusI(SgOps->SMx[IxGen[iGen]].s.R, &GenRmI[iRmI * 9], 0);
        iRmI++;
      }
    }
    else {
      range1(iGen, nGen) {
        if (CB_RMx(&GenRmI[iRmI * 9],
                   Z2PCBMx[0].s.R,
                   SgOps->SMx[IxGen[iGen]].s.R,
                   Z2PCBMx[1].s.R) != 0) return -1;
        range3(i, 0, 9, 4) GenRmI[iRmI * 9 + i] -= 1;
        iRmI++;
      }
    }

    if (iRmI * 3 != nrGenRmI) return IE(-1);

    return nrGenRmI;
  }


  static int GetContNullSpace(const T_SgOps *SgOps, int IxGen[2], int nGen,
                              T_ssVM ssVM[3])
  {
    int  GenRmI[3 * 3 * 3];
    int  nrGenRmI, RankGenRmI;
    int  nIndep, iIndep, IxIndep[3], Sol[4][3];
    int  n_ssVM;

        nrGenRmI = ConstructGenRmI(SgOps, NULL, IxGen, nGen, GenRmI);
    if (nrGenRmI < 0) return IE(-1);

        RankGenRmI = iRowEchelonFormT(GenRmI, nrGenRmI, 3, NULL, 0);
    if (RankGenRmI < 0 || RankGenRmI > 3)
      return IE(-1);

    n_ssVM = 3 - RankGenRmI;

        nIndep = iRESetIxIndep(GenRmI, RankGenRmI, 3, IxIndep, 3);
    if (nIndep < 0) return IE(-1);

    if (nIndep != 2)
    {
      range1(iIndep, nIndep)
      {
        ssVM[iIndep].V[IxIndep[iIndep]] = 1;

        if (iREBacksubst(GenRmI, NULL, RankGenRmI, 3,
                         ssVM[iIndep].V, NULL) < 1)
          return IE(-1);

        ssVM[iIndep].M = 0;
      }
    }
    else
    {
      if (SolveHomRE1(GenRmI, IxIndep, Sol) != 0) return -1;

      qsort((void *) Sol, 4, sizeof (*Sol),
            (int (*)(const void *, const void *)) CmpOLen2);

      range1(iIndep, 2) {
        MemCpy(ssVM[iIndep].V, Sol[iIndep], 3);
        ssVM[iIndep].M = 0;
      }
    }

    return n_ssVM;
  }


  static int nDLoopStep(int *i, int n, int Low, int High)
  {
    int  p, l;

    p = l = n - 1;

    for (; p >= 0;)
    {
          i[p]++;
      if (i[p] > High)
        p--;
      else if (p < l)
        i[++p] = Low - 1;
      else
        return 1;
    }

    return 0;
  }


  static void UpdateBestZ(int OrigZf[mDiscrGr][3], int nDiscrGr,
                          int BestZf[mDiscrGr][3], int BestM[mDiscrGr],
                          int BestZc[mDiscrGr][3],
                          int Shift[3], int LTBF)
  {
    int  iDG, c, i;
    int  Zf[3], M;
    int  Zc[3];

    range2(iDG, 1, nDiscrGr)
    {
      rangei(3) Zf[i] = iModPositive(Shift[i] + OrigZf[iDG][i], LTBF);
      MemCpy(Zc, Zf, 3);
      M = CancelBFGCD(Zc, 3, LTBF);

      rangei(3) {
        if (Zf[i]) {
          c = CmpOLen2(BestZc[iDG], Zc);
          if (c > 0 || (c == 0 && BestM[iDG] > M)) {
            MemCpy(BestZf[iDG], Zf, 3);
            MemCpy(BestZc[iDG], Zc, 3);
            BestM[iDG] = M;
          }
          break;
        }
      }
    }
  }


  static int BestVect(const T_SgOps *SgOps,
                      const T_ssVM *ssVM, int n_ssVM,
                      int DTBF, T_DiscrGr DiscrGr[mDiscrGr], int nDiscrGr)
  {
    int  iLTr, i_ssVM, iDG, gcd, i, j;
    int  fGrd, LTBF, f[2];
    int  OrigZf[mDiscrGr][3];
    int  BestZf[mDiscrGr][3], BestM[mDiscrGr];
    int  BestZc[mDiscrGr][3];
    int  LTr[3][3];

    fGrd = 1;
    LTBF = 1;

    range2(iDG, 1, nDiscrGr) {
      rangei(3) {
        gcd = iGCD(DiscrGr[iDG].Z[i], DTBF * CRBF);
        LTBF = iLCM(LTBF, (DTBF * CRBF) / gcd);
      }
    }

    range2(iLTr, 1, SgOps->nLTr) {
      rangei(3) {
        gcd = iGCD(SgOps->LTr[iLTr].v[i], STBF);
        LTBF = iLCM(LTBF, STBF / gcd);
      }
    }

    range1(i_ssVM, n_ssVM)
      rangei(3) fGrd = iLCM(fGrd, ssVM[i_ssVM].V[i]);

    LTBF *= fGrd;

    if (LTBF > 6) LTBF = iLCM(LTBF,  6);
    else          LTBF = iLCM(LTBF, 12);

    if (SgOps->nLTr == 1 && n_ssVM == 0) return 0;

    range2(iDG, 1, nDiscrGr) {
      if (ChangeBaseFactor(DiscrGr[iDG].Z, DTBF * CRBF,
                              OrigZf[iDG],        LTBF, 3) != 0)
        return IE(-1);
      rangei(3) OrigZf[iDG][i] = iModPositive(OrigZf[iDG][i], LTBF);
      MemCpy(BestZf[iDG], OrigZf[iDG], 3);
      MemCpy(BestZc[iDG], OrigZf[iDG], 3);
      BestM[iDG] = CancelBFGCD(BestZc[iDG], 3, LTBF);
    }

    if (n_ssVM > 2) return IE(-1);

    range1(iLTr, SgOps->nLTr)
    {
      if (ChangeBaseFactor(SgOps->LTr[iLTr].v, STBF, LTr[0], LTBF, 3) != 0)
        return IE(-1);

      rangei(n_ssVM) f[i] = 0;

      do
      {
        rangei(n_ssVM)
          range1(j, 3) LTr[i + 1][j] = LTr[i][j] + f[i] * ssVM[i].V[j];

        UpdateBestZ(OrigZf, nDiscrGr, BestZf, BestM, BestZc,
                    LTr[n_ssVM], LTBF);
      }
      while (nDLoopStep(f, n_ssVM, 0, LTBF - 1));
    }

    range2(iDG, 1, nDiscrGr)
      if (ChangeBaseFactor(BestZf[iDG],           LTBF,
                           DiscrGr[iDG].Z, DTBF * CRBF, 3) != 0)
        return IE(-1);

    return 0;
  }


  static int SelectDiscrete(int LTBF, int nDiscrGr, T_DiscrGr DiscrGr[mDiscrGr],
                            int mIx, int Ix[3])
  {
    int    nIx, iIx;
    T_LTr  LLTr[mDiscrGr];
    int    nLLTr;

    if (nDiscrGr == 1) return 0;

    for (nIx = 1; nIx <= nDiscrGr - 1 && nIx <= mIx; nIx++)
    {
      range1(iIx, nIx) Ix[iIx] = iIx;

      do
      {
        ResetLLTr(LLTr, &nLLTr);
        range1(iIx, nIx)
          if (ExpLLTr(LTBF, mDiscrGr,
                      LLTr, &nLLTr, DiscrGr[Ix[iIx] + 1].P) < 0)
            return IE(-1);

        if (nLLTr >  nDiscrGr) return IE(-1);
        if (nLLTr == nDiscrGr) return nIx;
      }
      while (NextOf_n_from_m(nDiscrGr - 1, nIx, Ix) != 0);
    }

    return IE(-1);
  }


  static int CmpDiscr(const T_DiscrGr *a, const T_DiscrGr *b)
  {
    return CmpiVect(a->Z, b->Z, 3);
  }


  static int Cmp_ssVM(const T_ssVM *a, const T_ssVM *b)
  {
    return CmpiVect(a->V, b->V, 3);
  }


  int Set_ss(const T_SgOps *SgOps, T_ssVM ssVM[3])
  {
    int        ir, i;
    int        nGen, IxGen[2];
    int        n_ssVM;
    T_RTMx     Z2PCBMx[2];
    int        nrSNF, DTBF, nd, id, d, f;
    int        SNF[3 * 3 * 3], Q[3 * 3], xp[3], x[3];
    int        nDiscrGr, iDG;
    T_LTr      LLTr[mDiscrGr];
    T_DiscrGr  DiscrGr[mDiscrGr];
    int        nIx, Ix[3];

    range1(ir, 3) rangei(3) ssVM[ir].V[i] = 0;
    range1(ir, 3) ssVM[ir].M = -1;

    // XXX new: P 1 -> ngen - 0
        nGen = SetAnyIxGen(SgOps, MGC_Undefined, IxGen);
    if (nGen < 0 || nGen > 2) return IE(-1);

        n_ssVM = GetContNullSpace(SgOps, IxGen, nGen, ssVM);
    if (n_ssVM < 0) return -1;
    if (n_ssVM == 3) return n_ssVM;

    if (GetZ2PCBMx(SgOps, Z2PCBMx) != 0) return -1;

        nrSNF = ConstructGenRmI(SgOps, Z2PCBMx, IxGen, nGen, SNF);
    if (nrSNF < 0) return IE(-1);

        nd = SmithNormalForm(SNF, nrSNF, 3, NULL, Q);
    if (nd < 0 || nd > 3) return IE(-1);

    DTBF = 1; range1(id, 3) DTBF = iLCM(DTBF, SNF[(nd + 1) * id]);

    ResetLLTr(LLTr, &nDiscrGr);

    range1(id, nd) {
      d = SNF[(nd + 1) * id];
      range2(f, 1, d) {
        rangei(3) xp[i] = 0;
        xp[id] = f * DTBF / d;
        iMxMultiply(x, xp, Q, 1, 3, 3);
        if (ExpLLTr(DTBF, mDiscrGr, LLTr, &nDiscrGr, x) < 0)
          return IE(-1);
      }
    }

    range1(iDG, nDiscrGr) {
      MemCpy(DiscrGr[iDG].P, LLTr[iDG].v, 3);
      RotMx_t_Vector(DiscrGr[iDG].Z, Z2PCBMx[1].s.R, DiscrGr[iDG].P, 0);
      rangei(3)
        DiscrGr[iDG].Z[i] = iModPositive(DiscrGr[iDG].Z[i], DTBF * CRBF);
    }

    if (BestVect(SgOps, ssVM, n_ssVM, DTBF, DiscrGr, nDiscrGr) != 0)
      return IE(-1);

    qsort((void *) DiscrGr, nDiscrGr, sizeof (*DiscrGr),
          (int (*)(const void *, const void *)) CmpDiscr);

        nIx = SelectDiscrete(DTBF, nDiscrGr, DiscrGr, 3 - n_ssVM, Ix);
    if (nIx < 0) return IE(-1);

    rangei(nIx) {
      if (n_ssVM >= 3) return IE(-1);
      MemCpy(ssVM[n_ssVM].V, DiscrGr[Ix[i] + 1].Z, 3);
      ssVM[n_ssVM].M = CancelBFGCD(ssVM[n_ssVM].V, 3, DTBF * CRBF);
      n_ssVM++;
    }

    qsort((void *) ssVM, n_ssVM, sizeof (*ssVM),
          (int (*)(const void *, const void *)) Cmp_ssVM);

    return n_ssVM;
  }
#endif

  StructureSeminvariant::StructureSeminvariant(const SgOps& sgo)
    : m_size(0)
  {
    for(std::size_t i = 0; i < 3; i++) m_VM[i].zero_out();

    // XXX new: P 1 -> ngen = 0
    detail::AnyGenerators Gen(sgo);

    int j;
    SgOps VfySgOps;
    for(j=0; j < sgo.nLTr(); j++) {
      VfySgOps.expandLTr(sgo.LTr(j));
    }
    if (Gen.ZInvT.isValid()) VfySgOps.expandInv(Gen.ZInvT);
    for(j=0; j < Gen.nGen; j++) {
      VfySgOps.expandSMx(Gen.ZGen[j]);
    }
    cctbx_assert(VfySgOps == sgo);
  }

} // namespace sgtbx
