// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Created: 24-Apr-2001 (R.W. Grosse-Kunstleve)
 */

#include <set>

#include <string>
#include <algorithm>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/reference.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  std::map<int, int> SgOps::CountRtypes() const
  {
    std::map<int, int> result;
    rangei(nSMx()) {
      result[m_SMx[i].Rpart().getRtype()]++;
    }
    return result;
  }

  tables::MatrixGroup::Code SgOps::getPointGroupType() const
  {
    using namespace tables::MatrixGroup;

    std::map<int, int> Counter = CountRtypes();

    if      (Counter[-3] + Counter[3] == 8) {
      if      (nSMx() == 12) {
        if (!isCentric()) return MGC_23;
        else              return MGC_m3b;
      }
      else if (nSMx() == 24) {
        if (!isCentric()) {
          if (Counter[ 4] == 6) return MGC_432;
          if (Counter[-4] == 6) return MGC_4b3m;
        }
        else return MGC_m3bm;
      }
    }
    else if (Counter[-6] + Counter[6] == 2) {
      if      (nSMx() ==  6) {
        if (!isCentric()) {
          if (Counter[ 6] == 2) return MGC_6;
          if (Counter[-6] == 2) return MGC_6b;
        }
        else return MGC_6_m;
      }
      else if (nSMx() == 12) {
        if (!isCentric()) {
          if (Counter[ 6] == 2) {
            if (Counter[ 2] == 7) return MGC_622;
            if (Counter[-2] == 6) return MGC_6mm;
          }
          else if (Counter[-6] == 2) return MGC_6bm2;
        }
        else return MGC_6_mmm;
      }
    }
    else if (Counter[-3] + Counter[3] == 2) {
      if      (nSMx() ==  3) {
        if (!isCentric()) return MGC_3;
        else return MGC_3b;
      }
      else if (nSMx() ==  6) {
        if (!isCentric()) {
          if (Counter[ 2] == 3) return MGC_32;
          if (Counter[-2] == 3) return MGC_3m;
        }
        else return MGC_3bm;
      }
    }
    else if (Counter[-4] + Counter[4] == 2) {
      if (nSMx() ==  4) {
        if (!isCentric()) {
          if (Counter[ 4] == 2) return MGC_4;
          if (Counter[-4] == 2) return MGC_4b;
        }
        else return MGC_4_m;
      }
      else if (nSMx() ==  8) {
        if (!isCentric()) {
          if (Counter[ 4] == 2) {
            if (Counter[ 2] == 5) return MGC_422;
            if (Counter[-2] == 4) return MGC_4mm;
          }
          else if (Counter[-4] == 2) return MGC_4bm2;
        }
        else return MGC_4_mmm;
      }
    }
    else if (Counter[-2] + Counter[2] == 3) {
      if (!isCentric()) {
        if (Counter[ 2] == 3) return MGC_222;
        if (Counter[-2] == 2) return MGC_mm2;
      }
      else return MGC_mmm;
    }
    else if (Counter[-2] + Counter[2] == 1) {
      if (!isCentric()) {
        if (Counter[ 2] == 1) return MGC_2;
        if (Counter[-2] == 1) return MGC_m;
      }
      else return MGC_2_m;
    }
    else if (nSMx() == 1) {
      if (!isCentric()) return MGC_1;
      else return MGC_1b;
    }

    throw cctbx_internal_error();
  }

  namespace ConstructCBOpRpart {

    tables::MatrixGroup::Code
    GetMatrixGroupType(const SgOps& StdSgOps,
                       tables::MatrixGroup::Code PointGroupType)
    {
      using namespace tables::MatrixGroup;

      int TwoFold = 0;
      int Mirror  = 0;

      if      (   PointGroupType == MGC_4bm2
               || PointGroupType == MGC_6bm2)
        TwoFold = 1;
      else if (StdSgOps.nLTr() == 1) {
        if      (PointGroupType == MGC_32)  TwoFold = 1;
        else if (PointGroupType == MGC_3m)  Mirror  = 1;
        else if (PointGroupType == MGC_3bm) TwoFold = Mirror = 1;
      }

      if (!(TwoFold || Mirror))
        return PointGroupType;

      for (int iSMx = 1; iSMx < StdSgOps.nSMx(); iSMx++) {
        int Rtype = StdSgOps[iSMx].Rpart().getRtype();
        cctbx_assert(Rtype != 0);
        if (   (Rtype ==  2 && TwoFold)
            || (Rtype == -2 && Mirror)) {
          RotMxInfo RI = StdSgOps[iSMx].Rpart().getInfo();
          const Vec3 EV_100 = {1, 0, 0};
          if (RI.EV() == EV_100) {
            if (PointGroupType == MGC_4bm2) return MGC_4b2m;
            if (PointGroupType == MGC_32)   return MGC_321;
            if (PointGroupType == MGC_3m)   return MGC_3m1;
            if (PointGroupType == MGC_3bm)  return MGC_3bm1;
            if (PointGroupType == MGC_6bm2) return MGC_6b2m;
          }
        }
      }

      if (PointGroupType == MGC_4bm2) return MGC_4bm2;
      if (PointGroupType == MGC_32)   return MGC_312;
      if (PointGroupType == MGC_3m)   return MGC_31m;
      if (PointGroupType == MGC_3bm)  return MGC_3b1m;
      if (PointGroupType == MGC_6bm2) return MGC_6bm2;

      throw cctbx_internal_error();
    }

    struct BasisBuilder {
      BasisBuilder(const tables::MatrixGroup::Code& LaueGroupType);
      void PickEigenVectors(const SgOps& RawSgOps);
      RotMx UniAxialBasis();
      RotMx TwoAxialBasis();
      const tables::CrystalSystem::Code CrystalSystem;
      int nWtd;
      int Ord[3];
      RotMx RMx[3];
      RotMxInfo RI[3];
    };

    BasisBuilder::BasisBuilder(const tables::MatrixGroup::Code& LaueGroupType)
      : CrystalSystem(LaueGroupType.CrystalSystem())
    {
      using namespace tables::MatrixGroup;

      if      (LaueGroupType == MGC_1b) {
        nWtd = 0;
      }
      else if (LaueGroupType == MGC_2_m) {
        nWtd = 1; Ord[0] = 2;
      }
      else if (LaueGroupType == MGC_mmm) {
        nWtd = 3; Ord[0] = 2; Ord[1] = 2; Ord[2] = 2;
      }
      else if (LaueGroupType == MGC_4_m) {
        nWtd = 1; Ord[0] = 4;
      }
      else if (LaueGroupType == MGC_4_mmm) {
        nWtd = 2; Ord[0] = 4; Ord[1] = 2;
      }
      else if (LaueGroupType == MGC_3b) {
        nWtd = 1; Ord[0] = 3;
      }
      else if (LaueGroupType == MGC_3bm) {
        nWtd = 2; Ord[0] = 3; Ord[1] = 2;
      }
      else if (LaueGroupType == MGC_6_m) {
        nWtd = 1; Ord[0] = 3;
      }
      else if (LaueGroupType == MGC_6_mmm) {
        nWtd = 2; Ord[0] = 3; Ord[1] = 2;
      }
      else if (LaueGroupType == MGC_m3b) {
        nWtd = 2; Ord[0] = 3; Ord[1] = 2;
      }
      else if (LaueGroupType == MGC_m3bm) {
        nWtd = 2; Ord[0] = 3; Ord[1] = 4;
      }
      else {
        throw cctbx_internal_error();
      }
    }

    void BasisBuilder::PickEigenVectors(const SgOps& RawSgOps)
    {
      int iWtd;
      for (iWtd = 0; iWtd < nWtd; iWtd++) RMx[iWtd] = RotMx(0);
      int nFound = 0;
      for (;;) {
        int Restart = 0;
        for (int iSMx = 0; iSMx < RawSgOps.nSMx(); iSMx++) {
          int RMxO = RawSgOps[iSMx].Rpart().getRtype();
          cctbx_assert(RMxO != 0);
          for (iWtd = 0; iWtd < nWtd; iWtd++) {
            if (   !RMx[iWtd].isValid()
                && (RMxO == Ord[iWtd] || -RMxO == Ord[iWtd])) {
              RI[iWtd] = RawSgOps[iSMx].Rpart().getInfo();
              if (RI[iWtd].SenseOfRotation() >= 0) {
                int UseThisSMx = 1;
                for (int jWtd = 0; jWtd < nWtd; jWtd++) {
                  if (RMx[jWtd].isValid() && RI[iWtd].EV() == RI[jWtd].EV()) {
                    if (std::abs(Ord[jWtd]) >= std::abs(RMxO)) {
                      UseThisSMx = 0;
                    }
                    else {
                      RMx[jWtd] = RotMx(0);
                      nFound--;
                      Restart = 1;
                    }
                    break;
                  }
                }
                if (UseThisSMx) {
                  Ord[iWtd] = RMxO;
                  RMx[iWtd] = RawSgOps[iSMx].Rpart();
                      nFound++;
                  if (nFound == nWtd)
                    return;
                }
              }
              break;
            }
          }
        }
        if (Restart == 0) break;
      }
      throw cctbx_internal_error();
    }

    void SolveHomRE1(const int REMx[3], const int IxIndep[2], Vec3 Sol[4])
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
        cctbx_assert(iREBacksubst(REMx, 0, 2, 3, Sol[iPV].elems, 0) > 0);
      }
    }

    RotMx BasisBuilder::UniAxialBasis()
    {
      if (Ord[0] < 0) {
        Ord[0] *= -1;
        RMx[0] *= -1;
        RI[0] = RMx[0].getInfo();
      }
      RotMx CumMx = RMx[0].accumulate(Ord[0]);
      cctbx_assert(iRowEchelonFormT(CumMx.elems, 3, 3, 0, 0) == 1);
      int IxIndep[2];
      cctbx_assert(iRESetIxIndep(CumMx.elems, 1, 3, IxIndep, 2) == 2);
      Vec3 Sol[4];
      SolveHomRE1(CumMx.elems, IxIndep, Sol);
      int nIx = 1; if (Ord[0] == 2) nIx++;
      int Ix[2];
      int i;
      for(i=0;i<nIx;i++) Ix[i] = i;
      RotMx TrialBasis;
      TrialBasis.setColumn(2, RI[0].EV());
      int MinDet = 0;
      RotMx result;
      do {
        for(i=0;i<nIx;i++) TrialBasis.setColumn(i, Sol[Ix[i]]);
        if (nIx == 1) TrialBasis.setColumn(1, RMx[0] * Sol[Ix[0]]);
        int Det = TrialBasis.det();
        if (Det != 0 && (MinDet == 0 || std::abs(MinDet) > std::abs(Det))) {
          MinDet = Det;
          result = TrialBasis;
        }
      }
      while (NextOf_n_from_m(4, nIx, Ix) != 0);
      cctbx_assert(MinDet != 0);
      if (MinDet < 0) result.swapColumns(0, 1);
      return result;
    }

    RotMx BasisBuilder::TwoAxialBasis()
    {
      RotMx result;
      result.setColumn(0, RI[1].EV());
      if (CrystalSystem == tables::CrystalSystem::Cubic) {
        for(int i=0;i<2;i++) {
          result.setColumn(i + 1, RMx[0] * result.getColumn(i));
        }
        if (result.det() < 0) result.swapColumns(0, 1);
      }
      else {
        result.setColumn(2, RI[0].EV());
        if (nWtd == 3) {
          result.setColumn(1, RI[2].EV());
          if (result.det() < 0) result.swapColumns(0, 1);
        }
        else {
          if (Ord[0] > 0) {
            result.setColumn(1,  RMx[0] * result.getColumn(0));
          }
          else {
            result.setColumn(1, -RMx[0] * result.getColumn(0));
          }
          cctbx_assert(result.det() > 0);
        }
      }
      return result;
    }

    RotMx StdBasis(const SgOps& RawSgOps, tables::MatrixGroup::Code MGC)
    {
      BasisBuilder BB(MGC.LaueGroupType());
      if (BB.nWtd == 0) return RotMx();
      BB.PickEigenVectors(RawSgOps);
      if (BB.nWtd == 1) return BB.UniAxialBasis();
      return BB.TwoAxialBasis();
    }

    bool m3b100Glide(const SgOps& WorkSgOps)
    {
      cctbx_assert(WorkSgOps.isCentric());
      for (int iSMx = 1; iSMx < WorkSgOps.nSMx(); iSMx++) {
        const RotMx& R = WorkSgOps[iSMx].Rpart();
        int Rtype = R.getRtype();
        cctbx_assert(Rtype != 0);
        if (std::abs(Rtype) == 2) {
          RotMxInfo RI = R.getInfo();
          const Vec3 EV_100 = {1, 0, 0};
          if (RI.EV() == EV_100) {
            int iInv = 0; if (Rtype == 2) iInv = 1;
            RTMx SMx = WorkSgOps(0, iInv, iSMx);
            TrVec wi = SMx.getIntrinsicPart();
            if (wi[2] % wi.BF() != 0) return true;
            return false;
          }
        }
      }
      throw cctbx_internal_error();
    }

    RotMx GetAdjRMx(const tables::MatrixGroup::Code& LG_MGC, char Z)
    {
      using namespace tables::MatrixGroup;

      if      (    LG_MGC == MGC_2_m) {
        // monoclinic unique c -> unique b
        return RotMx(  0,  1,  0,
                       0,  0,  1,
                       1,  0,  0,  1);
      }
      if ((LG_MGC == MGC_4_m || LG_MGC == MGC_4_mmm) && Z == 'C') {
        // C -> P
        return RotMx(  1,  1,  0,
                       1, -1,  0,
                       0,  0, -1,  1);
      }
      if ((LG_MGC == MGC_4_m || LG_MGC == MGC_4_mmm) && Z == 'F') {
        // F -> I
        return RotMx(  1,  1,  0,
                      -1,  1,  0,
                       0,  0,  1,  1);
      }
      if ((LG_MGC == MGC_3b || LG_MGC == MGC_3bm) && Z == 'Q') {
        // reverse -> obverse
        return RotMx( -1,  0,  0,
                       0, -1,  0,
                       0,  0,  1,  1);
      }
      if ((LG_MGC == MGC_3bm || LG_MGC == MGC_6_mmm) && Z == 'H') {
        // H -> P
        return RotMx(  1,  1,  0,
                      -1,  2,  0,
                       0,  0,  1,  1);
      }
      return RotMx(0);
    }

  } // namespace ConstructCBOpRpart
} // namespace sgtbx

namespace sgtbx {
  namespace ConstructCBOpTpart {

    struct Generators {
      Generators(const SgOps& WorkSgOps,
                 const tables::MatrixGroup::Code& PG_MGC);
      void setPrimitive();
      ChOfBasisOp Z2POp;
      TrVec ZInvT;
      TrVec PInvT;
      int nGen;
      RTMx ZGen[2];
      RTMx PGen[2];
    };

    Generators::Generators(const SgOps& WorkSgOps,
                           const tables::MatrixGroup::Code& PG_MGC)
      : nGen(0)
    {
      using namespace tables::CrystalSystem;

      const Vec3 EV_001 = { 0, 0, 1};
      const Vec3 EV_100 = { 1, 0, 0};
      const Vec3 EV_110 = { 1, 1, 0};
      const Vec3 EV_m10 = {-1, 1, 0};
      const Vec3 EV_111 = { 1, 1, 1};

      Z2POp = WorkSgOps.getZ2POp();

      ZInvT = WorkSgOps.InvT(true);
      PInvT = TrVec(0);
      int i;
      for(i=0;i<2;i++) ZGen[i] = RTMx(0, 0);
      for(i=0;i<2;i++) PGen[i] = RTMx(0, 0);

      int PrincipalProperOrder = 0;

      switch (PG_MGC.CrystalSystem())
      {
        case Triclinic:
          break;

        case Monoclinic:
          ZGen[0] = WorkSgOps[1];
          nGen = 1;
          break;

        case Orthorhombic:
          for(i = 1; i < WorkSgOps.nSMx(); i++) {
            RotMxInfo RI = WorkSgOps[i].Rpart().getInfo();
            if      (RI.EV() == EV_001) { ZGen[0] = WorkSgOps[i]; nGen++; }
            else if (RI.EV() == EV_100) { ZGen[1] = WorkSgOps[i]; nGen++; }
          }
          cctbx_assert(nGen == 2);
          break;

        case Tetragonal:
                                     PrincipalProperOrder = 4;
        case Trigonal:
          if (!PrincipalProperOrder) PrincipalProperOrder = 3;
        case Hexagonal:
          if (!PrincipalProperOrder) PrincipalProperOrder = 6;

          for(i = 1; i < WorkSgOps.nSMx(); i++) {
            RotMxInfo RI = WorkSgOps[i].Rpart().getInfo();
            if (std::abs(RI.Rtype()) == PrincipalProperOrder) {
              if (RI.SenseOfRotation() > 0) {
                ZGen[0] = WorkSgOps[i]; nGen++;
              }
            }
            else if (PrincipalProperOrder == 4) {
              if (RI.EV() == EV_100) { ZGen[1] = WorkSgOps[i]; nGen++; }
            }
            else if (PrincipalProperOrder == 3) {
              if      (RI.EV() == EV_m10) { ZGen[1] = WorkSgOps[i]; nGen++; }
              else if (RI.EV() == EV_110) { ZGen[1] = WorkSgOps[i]; nGen++; }
            }
            else { // PrinicipalProperOrder == 6
              if (RI.EV() == EV_m10) { ZGen[1] = WorkSgOps[i]; nGen++; }
            }
          }
          cctbx_assert(nGen == 1 || nGen == 2);
          for (i=0;i<nGen;i++) cctbx_assert(ZGen[i].isValid());
          break;

        case Cubic:
          for(i = 1; i < WorkSgOps.nSMx(); i++) {
            RotMxInfo RI = WorkSgOps[i].Rpart().getInfo();
            if      (std::abs(RI.Rtype()) == 4) {
              if (RI.SenseOfRotation() > 0 && RI.EV() == EV_001) {
                if (!ZGen[0].isValid()) nGen++;
                ZGen[0] = WorkSgOps[i];
              }
            }
            else if (std::abs(RI.Rtype()) == 2) {
              if (!ZGen[0].isValid() && RI.EV() == EV_001) {
                ZGen[0] = WorkSgOps[i]; nGen++;
              }
            }
            else if (std::abs(RI.Rtype()) == 3) {
              if (RI.SenseOfRotation() > 0 && RI.EV() == EV_111) {
                cctbx_assert(!ZGen[1].isValid());
                ZGen[1] = WorkSgOps[i]; nGen++;
              }
            }
          }
          cctbx_assert(nGen == 1 || nGen == 2);
          for (i=0;i<nGen;i++) cctbx_assert(ZGen[i].isValid());
          break;

        default:
          throw cctbx_internal_error();
      }

      // Tidy generators
      if (ZInvT.isValid()) {
        for (i = 0; i < nGen; i++) {
          if (ZGen[i].Rpart().det() < 0) {
            ZGen[i] = ZGen[i].pre_multiply_InvT(ZInvT);
          }
        }
      }
      for (i = 0; i < nGen; i++) ZGen[i].ModPositive();
    }

    void Generators::setPrimitive()
    {
      for (int i = 0; i < nGen; i++) {
        PGen[i] = Z2POp(ZGen[i]);
        PGen[i].ModPositive();
      }
      if (ZInvT.isValid()) {
        PInvT = Z2POp(ZInvT, -1);
        PInvT.ModPositive();
      }
    }

    bool SolveInhomModZ(int *M, int nr, int nc, int *b, int BF, int *x)
    {
      const int maxr = 9;
      const int maxc = 6;

      cctbx_assert(nr <= maxr);
      cctbx_assert(nc <= maxc);

      int P[maxr * maxr], Q[maxc * maxc];
      int nd = SmithNormalForm(M, nr, nc, P, Q);
      cctbx_assert(nd >= 0 && nd <= nc);

      int Pb[maxr];
      MatrixLite::multiply<int>(P, b, nr, nr, 1, Pb);
      int i;
      for(i = nd; i < nr; i++) {
        if (Pb[i] % BF != 0) return false;
      }

      if (x) {
        int xp[maxc];
        for(i=0;i<nc;i++) {
          xp[i] = 0;
          int d = M[i * nd + i];
          if (d) {
            cctbx_assert(Pb[i] % d == 0);
            xp[i] = Pb[i] / d;
          }
        }
        MatrixLite::multiply<int>(xp, Q, 1, nc, nc, x);
      }

      return true;
    }

    TrVec FindOriginShift(const Generators& TabGenerators,
                          const Generators& TstGenerators,
                          int TBF)
    {
      /*    (I|K)(R|T)(I|-K)=(R|S)
         => S=-RK+T+K=-(R-I)K+T
         => S=-(R-I)K+T
         => (R-I)K=T-S
         => (R-I)^-1(T-S)=K
       */

      RotMx RmI[3];
      TrVec DeltaT[3];
      int i;
      for (i = 0; i < TabGenerators.nGen; i++) {
        cctbx_assert(   TstGenerators.PGen[i].Rpart()
                     == TabGenerators.PGen[i].Rpart());
        RmI[i] = TstGenerators.PGen[i].Rpart() - RotMx();
        DeltaT[i] = (  TstGenerators.PGen[i].Tpart()
                     - TabGenerators.PGen[i].Tpart()).newBaseFactor(TBF);
      }
      cctbx_assert(   TstGenerators.PInvT.isValid()
                   == TabGenerators.PInvT.isValid());
      if (TstGenerators.PInvT.isValid()) {
        RmI[i] = RotMx(1, -2);
        DeltaT[i] = (  TstGenerators.PInvT
                     - TabGenerators.PInvT).newBaseFactor(TBF);
        i++;
      }
      const int nGen = i;

      const int nrSNF = nGen * 3;
      int SNF[9 * 3], V[3 * 3];
      for (int iGen = 0; iGen < nGen; iGen++) {
        for(i=0;i<9;i++) SNF[iGen * 9 + i] = RmI[iGen][i];
        for(i=0;i<3;i++) V[iGen * 3 + i] = DeltaT[iGen][i];
      }

      Vec3 x;
      if (SolveInhomModZ(SNF, nrSNF, 3, V, TBF, x.elems)) {
        TrVec CBT = TstGenerators.Z2POp.InvM().Rpart() * TrVec(x, TBF);
        return CBT.newBaseFactor(TBF);
      }

      return TrVec(0);
    }

    ChOfBasisOp MatchGenerators(int RBF, int TBF,
                                const SgOps& WorkSgOps,
                                const tables::MatrixGroup::Code& PG_MGC,
                                char TabZ,
                                const Generators& TabGenerators)
    {
      if (TabGenerators.nGen == 0 && !TabGenerators.ZInvT.isValid())
        return ChOfBasisOp(RBF, TBF); // space group P 1

      RotMx RMx_2fold(0);
      RotMx RMx_3fold(0);

      if      (PG_MGC.CrystalSystem() == tables::CrystalSystem::Monoclinic) {
        RMx_2fold = RotMx( // 2 [101]
          0,  0,  1,
          0, -1,  0,
          1,  0,  0,  1);
        RMx_3fold = RotMx( // 3 [010]
         -1,  0,  1,
          0,  1,  0,
         -1,  0,  0,  1);
      }
      else if (PG_MGC.CrystalSystem() == tables::CrystalSystem::Orthorhombic) {
        RMx_2fold = RotMx( // 2 [110]
          0,  1,  0,
          1,  0,  0,
          0,  0, -1,  1);
        RMx_3fold = RotMx( // 3 [111]
          0,  0,  1,
          1,  0,  0,
          0,  1,  0,  1);
      }

      if (RMx_2fold.isValid())
      {
        const ChOfBasisOp
        CBOp_2fold(RTMx(RMx_2fold.newBaseFactor(RBF), TBF));
        const ChOfBasisOp
        CBOp_3fold(RTMx(RMx_3fold.newBaseFactor(RBF), TBF));

        ChOfBasisOp CBOp(RBF, TBF);
        for (int i2fold = 0; i2fold < 2; i2fold++) {
          if (i2fold) CBOp = CBOp_2fold;
          for (int i3fold = 0; i3fold < 3; i3fold++) {
            if (i3fold) CBOp.update(CBOp_3fold);
            SgOps TstSgOps = WorkSgOps.ChangeBasis(CBOp);
            const char TstZ = TstSgOps.getConventionalCentringTypeSymbol();
            cctbx_assert(TstZ != '\0' && TstZ != 'Q');
            if (TstZ != TabZ)
              continue;
            Generators TstGenerators(TstSgOps, PG_MGC);
            cctbx_assert(TstGenerators.nGen == TabGenerators.nGen);
            if (    TabGenerators.nGen != 2
                ||  (      TabGenerators.ZGen[0].Rpart()[8]
                        == TstGenerators.ZGen[0].Rpart()[8]
                     &&    TabGenerators.ZGen[1].Rpart()[0]
                        == TstGenerators.ZGen[1].Rpart()[0])) {
              TstGenerators.setPrimitive();
              TrVec CBT = FindOriginShift(TabGenerators, TstGenerators, TBF);
              if (CBT.isValid()) {
                CBOp.update(CBT);
                return CBOp;
              }
            }
          }
        }
      }
      else {
        Generators TstGenerators(WorkSgOps, PG_MGC);
        cctbx_assert(TstGenerators.nGen == TabGenerators.nGen);
        TstGenerators.setPrimitive();
        TrVec CBT = FindOriginShift(TabGenerators, TstGenerators, TBF);
        if (CBT.isValid()) {
          return ChOfBasisOp(RTMx(CBT, RBF));
        }
      }

      return ChOfBasisOp(0, 0);
    }

    void TidyCBOpT(const Generators& TargetGenerators,
                   const SgOps& GivenSgOps,
                   const tables::MatrixGroup::Code& PG_MGC,
                   ChOfBasisOp& TrialCBOp)
    {
      // set translation parts to zero
      TrialCBOp = ChOfBasisOp(RTMx(TrialCBOp.M().Rpart(),
                                   TrialCBOp.M().Tpart().Null()),
                              RTMx(TrialCBOp.InvM().Rpart(),
                                   TrialCBOp.InvM().Tpart().Null()));

      // done if space group is P 1
      if (GivenSgOps.nSMx() == 1 && !GivenSgOps.isCentric()) return;

      SgOps TransformedSgOps = GivenSgOps.ChangeBasis(TrialCBOp);
      Generators TransformedGenerators(TransformedSgOps, PG_MGC);
      TransformedGenerators.setPrimitive();
      TrVec
      CBT = FindOriginShift(TargetGenerators, TransformedGenerators,
                            TrialCBOp.M().TBF());
      cctbx_assert(CBT.isValid());
      TrialCBOp.update(CBT.ModShort());
    }

    class CmpChOfBasisMx {
      public:
        CmpChOfBasisMx() : CmpR(9), CmpT(3) {}
        bool operator()(const RTMx& a, const RTMx& b) const
        {
          const RotMx& aR = a.Rpart(); const TrVec& aT = a.Tpart();
          const RotMx& bR = b.Rpart(); const TrVec& bT = b.Tpart();

          bool ba = aR.isUnit();
          bool bb = bR.isUnit();
          if ( ba && !bb) return true;
          if (!ba &&  bb) return false;

          ba = aT.isNull();
          bb = bT.isNull();
          if ( ba && !bb) return true;
          if (!ba &&  bb) return false;

          int i;
          int na = 0; for(i=0;i<9;i++) if (aR[i] == 0) na++;
          int nb = 0; for(i=0;i<9;i++) if (bR[i] == 0) nb++;
          if (na > nb) return true;
          if (na < nb) return false;

          na = 0; for(i=0;i<9;i++) if (std::abs(aR[i]) == aR.BF()) na++;
          nb = 0; for(i=0;i<9;i++) if (std::abs(bR[i]) == bR.BF()) nb++;
          if (na > nb) return true;
          if (na < nb) return false;

          na = 0; for(i=0;i<9;i++) if (aR[i] > 0) na++;
          nb = 0; for(i=0;i<9;i++) if (bR[i] > 0) nb++;
          if (na > nb) return true;
          if (na < nb) return false;

          if (CmpT(aT.elems, bT.elems)) return true;
          if (CmpT(bT.elems, aT.elems)) return false;

          return CmpR(bR.elems, aR.elems);
        }
      private:
        const CmpiVect CmpR;
        const CmpiVect CmpT;
    };

    /* Use affine normalizer operations to find the "best"
       change-of-basis matrix.
       For a given space group representation the resulting
       change-of-basis matrix should always be identical,
       independent of the order of symmetry operations
       and indepenent of the RawCBOp.
     */
    ChOfBasisOp FindBestCBOp(const SgOps& GivenSgOps,
                             const tables::MatrixGroup::Code& PG_MGC,
                             int SgNumber,
                             const ChOfBasisOp& RefCBOp,
                             const SgOps& TargetSgOps,
                             const Generators& TargetGenerators,
                             const ChOfBasisOp& RawCBOp)
    {
      std::vector<RTMx>
      AddlG = ReferenceSettings::GetNormAddlG(SgNumber, true, true, true);
      ChOfBasisOp CBOp = RefCBOp.swap();
      SgOps NormSgOps = TargetSgOps;
      int OldOrderP = NormSgOps.OrderP();
      for(int i=0;i<AddlG.size();i++) {
        RTMx M = CBOp(AddlG[i]);
        NormSgOps.expandSMx(M);
        int NewOrderP = NormSgOps.OrderP();
        cctbx_assert(   OldOrderP * OrderOfRtype(M.Rpart().getRtype())
                     == NewOrderP);
        OldOrderP = NewOrderP;
      }

      ChOfBasisOp BestCBOp = RawCBOp;
      bool CmpInv = false;
      if (BestCBOp.M().Rpart().det() < BestCBOp.InvM().Rpart().det())
        CmpInv = true;
      CmpChOfBasisMx CmpCBMx;

      for (int iInv = 0; iInv < NormSgOps.fInv(); iInv++)
      for (int iSMx = 0; iSMx < NormSgOps.nSMx(); iSMx++)
      {
        RTMx M = NormSgOps(0, iInv, iSMx);
        if (M.Rpart().det() < 0) continue;
        ChOfBasisOp NormCBOp(M.newBaseFactors(RawCBOp.M()));
        ChOfBasisOp TrialCBOp = NormCBOp * RawCBOp;
        if (SgNumber < 3 || SgNumber > 15) {
          TidyCBOpT(TargetGenerators, GivenSgOps, PG_MGC, TrialCBOp);
          if (!CmpCBMx(BestCBOp.select(CmpInv), TrialCBOp.select(CmpInv)))
            BestCBOp = TrialCBOp;
        }
        else {
          const int RBF = RawCBOp.M().RBF();
          const int TBF = RawCBOp.M().TBF();
          int r00, r22;
          ReferenceSettings::GetMonoAffNormTrialRanges(
            TrialCBOp.M().Rpart(), r00, r22);
          Mx33 M, InvM;
          M.assign(0);
          InvM.assign(0);
#define loop(i, r) \
          for (M[i] = -r * RBF; M[i] <= r * RBF; M[i] += RBF)
          loop(0, r00)
          loop(2, r22)
          loop(6, r00)
          loop(8, r22) {
#undef loop
            if (!ReferenceSettings::CheckMonoAffNormRestrictions(
                                         SgNumber, RotMx(M, RBF)))
              continue;
            M[4] = M[0] * M[8] - M[2] * M[6];
            if (   M[4] != -RBF * RBF
                && M[4] !=  RBF * RBF) continue;
            M[4] /= RBF;
            int f = M[4] / RBF;
            InvM[0] =  f * M[8];
            InvM[2] = -f * M[2];
            InvM[4] =      M[4];
            InvM[6] = -f * M[6];
            InvM[8] =  f * M[0];
            ChOfBasisOp M_TrialCBOp = ChOfBasisOp(
                RTMx(RotMx(M, RBF), TBF), RTMx(RotMx(InvM, RBF), TBF))
              * TrialCBOp;
            TidyCBOpT(TargetGenerators, GivenSgOps, PG_MGC, M_TrialCBOp);
            if (!CmpCBMx(BestCBOp.select(CmpInv),M_TrialCBOp.select(CmpInv)))
              BestCBOp = M_TrialCBOp;
          }
        }
      }
      return BestCBOp;
    }

  } // namespace ConstructCBOpTpart

  SpaceGroupType SgOps::getSpaceGroupType(bool TidyCBOp,
                                          int RBF, int TBF) const
  {
    tables::MatrixGroup::Code PG_MGC = getPointGroupType();
    tables::MatrixGroup::Code LG_MGC = PG_MGC.LaueGroupType();
    tables::CrystalSystem::Code CrystalSystem = PG_MGC.CrystalSystem();

    SgOps ZPointGroup = toZPointGroup(); // set translation parts to zero
    SgOps WorkSgOps = ZPointGroup;
    ChOfBasisOp TotCBOp(RBF, TBF);
    char Z = '\0';
    int RunAwayCounter = 0;

    do // align the point group operations
    {
      cctbx_assert(RunAwayCounter++ < 10);

      ChOfBasisOp AddCBOp = WorkSgOps.getZ2POp();
      TotCBOp.update(AddCBOp.newBaseFactors(RBF, TBF));
      WorkSgOps = ZPointGroup.ChangeBasis(TotCBOp);
      cctbx_assert(WorkSgOps.nLTr() == 1);

      RotMx Basis = ConstructCBOpRpart::StdBasis(WorkSgOps, LG_MGC);
      Basis = Basis.newBaseFactor(RBF);
      AddCBOp = ChOfBasisOp(RTMx(Basis, TBF)).swap();
      TotCBOp.update(AddCBOp);
      WorkSgOps = ZPointGroup.ChangeBasis(TotCBOp);
      Z = WorkSgOps.getConventionalCentringTypeSymbol();

      RotMx AdjRMx = ConstructCBOpRpart::GetAdjRMx(LG_MGC, Z);
      if (AdjRMx.isValid()) {
        AdjRMx = AdjRMx.newBaseFactor(RBF);
        AddCBOp = ChOfBasisOp(RTMx(AdjRMx, TBF));
        TotCBOp.update(AddCBOp);
        WorkSgOps = ZPointGroup.ChangeBasis(TotCBOp);
        Z = WorkSgOps.getConventionalCentringTypeSymbol();
      }
    }
    while (Z == '\0');
    cctbx_assert(Z != 'Q');

    WorkSgOps = ChangeBasis(TotCBOp); // transform original symmetry operations

    if (   PG_MGC == tables::MatrixGroup::MGC_m3b
        && Z == 'P'
        && ConstructCBOpRpart::m3b100Glide(WorkSgOps)) {
      // rotate by 90 degree (4 [0 0 1])
      RotMx R_4_001(  0, -1,  0,
                      1,  0,  0,
                      0,  0,  1,  1);
      ChOfBasisOp AddCBOp(RTMx(R_4_001.newBaseFactor(RBF), TBF));
      TotCBOp.update(AddCBOp);
      WorkSgOps = ChangeBasis(TotCBOp);
    }

    tables::MatrixGroup::Code
    MGC = ConstructCBOpRpart::GetMatrixGroupType(WorkSgOps, PG_MGC);
    cctbx_assert(PG_MGC == MGC.PointGroupType());

    bool MatchSymCType = (
             CrystalSystem != tables::CrystalSystem::Monoclinic
      && (   CrystalSystem != tables::CrystalSystem::Orthorhombic
          || (Z == 'I' || Z == 'F')));

    for (int SgNumber = 1; SgNumber <= 230; SgNumber++)
    {
      const char*
      HallSymbol = tables::ReferenceSettings::HallSymbols[SgNumber];

      if (MatchSymCType && Z != HallSymbol[1])
        continue;

      if ((Z == 'P') != (HallSymbol[1] == 'P'))
        continue;

      if (tables::ReferenceSettings::MatrixGroupCodes[SgNumber] != MGC)
        continue;

      SgOps TabSgOps;
      try {
        TabSgOps = SgOps(HallSymbol, true);
      }
      catch (const error&) {
        throw cctbx_internal_error();
      }

      if (TabSgOps.nLTr() != WorkSgOps.nLTr())
        continue;

      ConstructCBOpTpart::Generators TabGenerators(TabSgOps, PG_MGC);
      TabGenerators.setPrimitive();
      ChOfBasisOp AddCBOp = ConstructCBOpTpart::MatchGenerators(
        RBF, TBF, WorkSgOps, PG_MGC, HallSymbol[1], TabGenerators);
      if (AddCBOp.isValid()) {
        TotCBOp.update(AddCBOp);
        if (TidyCBOp) {
          TotCBOp = ConstructCBOpTpart::FindBestCBOp(
            *this, PG_MGC,
            SgNumber, ChOfBasisOp(RBF, TBF), TabSgOps, TabGenerators, TotCBOp);
        }
        return SpaceGroupType(SgNumber, TotCBOp);
      }
    }
    throw cctbx_internal_error();
  }

  std::string SgOps::BuildHallSymbol(const SpaceGroupType& SgType,
                                     bool TidyCBOp) const
  {
    cctbx_assert(0 < SgType.SgNumber() && SgType.SgNumber() <= 230);
    std::string
    HallSymbol(tables::ReferenceSettings::HallSymbols[SgType.SgNumber()]);
    parse_string ps(HallSymbol);
    SgOps TargetSgOps;
    ChOfBasisOp RefCBOp;
    TargetSgOps.ParseHallSymbolCBOp(ps, RefCBOp, true);
    ChOfBasisOp CBOp = RefCBOp;
    if (CBOp.isValid()) {
      CBOp = CBOp.newBaseFactors(SgType.CBOp());
      CBOp = CBOp.swap() * SgType.CBOp();
    }
    else {
      CBOp = SgType.CBOp();
    }
    if (TidyCBOp) {
      tables::MatrixGroup::Code PG_MGC = TargetSgOps.getPointGroupType();
      ConstructCBOpTpart::Generators TargetGenerators(TargetSgOps, PG_MGC);
      TargetGenerators.setPrimitive();
      CBOp = ConstructCBOpTpart::FindBestCBOp(
        *this, PG_MGC,
        SgType.SgNumber(), RefCBOp, TargetSgOps, TargetGenerators, CBOp);
    }
    else {
      CBOp.ModShort();
    }

    std::string::size_type par = HallSymbol.find(" (");
    if (par != std::string::npos) {
      HallSymbol.resize(par);
    }

    if (!CBOp.isUnit()) {
      HallSymbol += " (" + CBOp.InvM().as_xyz() + ")";
    }

    return HallSymbol;
  }

  std::string SgOps::BuildHallSymbol(bool TidyCBOp) const {
    return BuildHallSymbol(getSpaceGroupType(), TidyCBOp);
  }

} // namespace sgtbx
