// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 02: start port of sglite/sgss.c (R.W. Grosse-Kunstleve)
 */

#include <algorithm>
#include <cctbx/array_family/loops.h>
#include <cctbx/math/loops.h>
#include <cctbx/sgtbx/select_generators.h>
#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/reference.h>
#include <cctbx/array_family/tiny_algebra.h>

namespace cctbx { namespace sgtbx {

  namespace seminvariant_cpp { // anonymous namespace causes
                               // warnings with some compilers

    inline void copy(const RotMx& source, int* target, std::size_t n) {
      for(std::size_t i=0;i<n;i++) target[i] = source[i];
    }

    af::tiny<int, 3 * 3*3>
    ConstructGenRmI(const detail::AnyGenerators& Gen, bool Primitive)
    {
      af::tiny<int, 3 * 3*3> result;
      for(std::size_t i=0;i<Gen.nGen;i++) {
        const RotMx R = Gen.ZGen[i].Rpart();
        if (!Primitive) {
          copy(R.minusUnit(), &result[i * 3*3], 3*3);
        }
        else {
          copy(Gen.Z2POp(R).minusUnit(), &result[i * 3*3], 3*3);
        }
      }
      if (Gen.ZInvT.isValid()) {
        copy(RotMx(1, -2), &result[Gen.nGen * 3*3], 3*3);
      }
      return result;
    }

    class CmpOLen2 {
      public:
        CmpOLen2() : CmpT(3) {}
        bool operator()(const af::int3& a, const af::int3& b) const
        {
          int i;
          int OLen2a = 0; for(i=0;i<3;i++) OLen2a += a[i] * a[i];
          int OLen2b = 0; for(i=0;i<3;i++) OLen2b += b[i] * b[i];
          if (OLen2a < OLen2b) return true;
          if (OLen2a > OLen2b) return false;
          return CmpT(a.begin(), b.begin());
        }
      private:
        const CmpiVect CmpT;
    };

    af::small<ssVM, 3>
    GetContNullSpace(const detail::AnyGenerators& Gen)
    {
      af::small<ssVM, 3> result;
      af::tiny<int, 3 * 3 * 3> GenRmI = ConstructGenRmI(Gen, false);
      int RankGenRmI = iRowEchelonFormT(
        GenRmI.begin(), Gen.nAll() * 3, 3, 0, 0);
      cctbx_assert(RankGenRmI >= 0 && RankGenRmI <= 3);
      int IxIndep[3];
      int nIndep = iRESetIxIndep(GenRmI.begin(), RankGenRmI, 3, IxIndep, 3);
      cctbx_assert(nIndep >= 0);
      if (nIndep != 2) {
        for (int iIndep = 0; iIndep < nIndep; iIndep++) {
          ssVM VM;
          if (Gen.nAll() == 0) VM.V.fill(0);
          VM.V[IxIndep[iIndep]] = 1;
          cctbx_assert(iREBacksubst(GenRmI.begin(), 0, RankGenRmI, 3,
                                    VM.V.begin(), 0) > 0);
          VM.M = 0;
          result.push_back(VM);
        }
      }
      else {
        af::tiny<af::int3, 4> Sol;
        SolveHomRE1(GenRmI.begin(), IxIndep, Sol.begin());
        std::sort(Sol.begin(), Sol.end(), CmpOLen2());
        for (int iIndep = 0; iIndep < 2; iIndep++) {
          ssVM VM;
          VM.V = Sol[iIndep];
          VM.M = 0;
          result.push_back(VM);
        }
      }
      return result;
    }

    void UpdateBestZ(const af::small<TrVec, 8>& OrigZf,
                     const TrVec& Shift,
                     af::small<TrVec, 8>& BestZf,
                     af::small<TrVec, 8>& BestZc)
    {
      for (int iDL = 1; iDL < OrigZf.size(); iDL++) {
        TrVec Zf = (OrigZf[iDL] + Shift).modPositive();
        TrVec Zc = Zf.cancel();
        for(int i=0;i<3;i++) {
          if (Zf[i]) {
            if (   CmpOLen2()(Zc.vec(), BestZc[iDL].vec())
                || (Zc.vec() == BestZc[iDL].vec()
                    && Zc.BF() < BestZc[iDL].BF())) {
              BestZf[iDL] = Zf;
              BestZc[iDL] = Zc;
            }
            break;
          }
        }
      }
    }

    void BestVectors(const SpaceGroup& SgOps,
                     const af::small<ssVM, 3>& ContinuousVM,
                     af::small<TrVec, 8>& DiscrZ)
    {
      if (SgOps.nLTr() == 1 && ContinuousVM.size() == 0) return;
      int LTBF = 1;
      int iDL;
      for (iDL = 1; iDL < DiscrZ.size(); iDL++) {
        int BF = DiscrZ[iDL].BF();
        for(int i=0;i<3;i++) {
          int g = gcd(DiscrZ[iDL][i], BF);
          LTBF = lcm(LTBF, BF / g);
        }
      }
      int iLTr;
      for (iLTr = 1; iLTr < SgOps.nLTr(); iLTr++) {
        int BF = SgOps.LTr(iLTr).BF();
        for(int i=0;i<3;i++) {
          int g = gcd(SgOps.LTr(iLTr)[i], BF);
          LTBF = lcm(LTBF, BF / g);
        }
      }
      int fGrd = 1;
      for (int iVM = 0; iVM < ContinuousVM.size(); iVM++) {
        for(int i=0;i<3;i++) {
          if (ContinuousVM[iVM].V[i]) fGrd = lcm(fGrd, ContinuousVM[iVM].V[i]);
        }
      }
      LTBF *= fGrd;
      cctbx_assert(LTBF > 0);
      if (LTBF > 6) LTBF = lcm(LTBF,  6);
      else          LTBF = lcm(LTBF, 12);
      af::small<TrVec, 8> OrigZf;
      af::small<TrVec, 8> BestZf;
      af::small<TrVec, 8> BestZc;
      for (iDL = 0; iDL < DiscrZ.size(); iDL++) {
        OrigZf.push_back(DiscrZ[iDL].newBaseFactor(LTBF));
        BestZf.push_back(OrigZf[iDL]);
        BestZc.push_back(BestZf[iDL].cancel());
      }
      cctbx_assert(ContinuousVM.size() < 3);
      af::small<int, 2> loop_end(ContinuousVM.size());
      std::fill(loop_end.begin(), loop_end.end(), LTBF);
      TrVec LTr[3];
      for (iLTr = 0; iLTr < SgOps.nLTr(); iLTr++) {
        LTr[0] = SgOps.LTr(iLTr).newBaseFactor(LTBF);
        af::nested_loop<af::small<int, 2> > loop(loop_end);
        do {
          for (int iVM = 0; iVM < ContinuousVM.size(); iVM++) {
            const af::small<int, 2>& f = loop();
            LTr[iVM + 1] = LTr[iVM] + f[iVM] * TrVec(ContinuousVM[iVM].V, LTBF);
          }
          UpdateBestZ(OrigZf, LTr[ContinuousVM.size()], BestZf, BestZc);
          loop.incr();
        }
        while (loop.over() == 0);
      }
      for (iDL = 1; iDL < DiscrZ.size(); iDL++) {
        DiscrZ[iDL] = BestZf[iDL].newBaseFactor(DiscrZ[iDL].BF());
      }
    }

    af::small<TrVec, 3>
    SelectDiscreteGenerators(const af::small<TrVec, 8>& DiscrP,
                             const af::small<TrVec, 8>& DiscrZ)
    {
      if (DiscrP.size() == 1) return af::small<TrVec, 3>();
      for (int nDVM = 1; nDVM <= DiscrP.size() - 1 && nDVM <= 3; nDVM++) {
        for (loop_n_from_m<3>
             Ix(DiscrP.size() - 1, nDVM); !Ix.over(); Ix.incr()) {
          TrOps DiscrGrp(DiscrP[0].BF());
          for (std::size_t iIx = 0; iIx < Ix.n(); iIx++) {
            DiscrGrp.expand(DiscrP[Ix[iIx] + 1]);
          }
          if (DiscrGrp.nVects() == DiscrP.size()) {
            af::small<TrVec, 3> result;
            for (std::size_t iIx = 0; iIx < Ix.n(); iIx++) {
              result.push_back(DiscrZ[Ix[iIx] + 1]);
            }
            return result;
          }
          cctbx_assert(DiscrGrp.nVects() < DiscrP.size());
        }
      }
      throw cctbx_internal_error();
    }

    class CmpDiscrZ {
      public:
        CmpDiscrZ() : CmpT(3) {}
        bool operator()(const TrVec& a, const TrVec& b) const {
          return CmpT(a.vec().begin(), b.vec().begin());
        }
      private:
        const CmpiVect CmpT;
    };

    class Cmp_ssVM {
      public:
        Cmp_ssVM() : CmpT(3) {}
        bool operator()(const ssVM& a, const ssVM& b) const {
          return CmpT(a.V.begin(), b.V.begin());
        }
      private:
        const CmpiVect CmpT;
    };

  } // namespace seminvariant_cpp

  using namespace seminvariant_cpp;

  StructureSeminvariant::StructureSeminvariant(const SpaceGroup& SgOps)
  {
    detail::AnyGenerators Gen(SgOps);
    m_VM = GetContNullSpace(Gen);
    if (m_VM.size() == 3) return; // space group P1
    af::tiny<int, 3 * 3 * 3> SNF = ConstructGenRmI(Gen, true);
    int Q[3 * 3];
    int nd = SmithNormalForm(SNF.begin(), Gen.nAll() * 3, 3, 0, Q);
    cctbx_assert(nd >=0 && nd <= 3);
    int id;
    int DTBF = 1;
    for (id = 0; id < nd; id++) DTBF = lcm(DTBF, SNF[(nd + 1) * id]);
    TrOps DiscrGrpP(DTBF);
    for (id = 0; id < nd; id++) {
      int d = SNF[(nd + 1) * id];
      for (int f = 1; f < d; f++) {
        af::int3 xp;
        xp.fill(0);
        xp[id] = f * DTBF / d;
        TrVec x(DTBF);
        MatrixLite::multiply<int>(Q, xp.begin(), 3, 3, 1, x.vec().begin());
        DiscrGrpP.expand(x);
      }
    }
    int iDL;
    af::small<TrVec, 8> DiscrZ;
    for (iDL = 0; iDL < DiscrGrpP.nVects(); iDL++) {
      DiscrZ.push_back(
        (Gen.Z2POp.InvM().Rpart() * DiscrGrpP[iDL]).modPositive());
    }
    BestVectors(SgOps, m_VM, DiscrZ);
    std::sort(DiscrZ.begin(), DiscrZ.end(), CmpDiscrZ());
    af::small<TrVec, 8> DiscrP;
    for (iDL = 0; iDL < DiscrZ.size(); iDL++) {
      DiscrP.push_back(
        (Gen.Z2POp.M().Rpart() * DiscrZ[iDL]).newBaseFactor(DTBF));
    }
    af::small<TrVec, 3>
    DiscrGen = SelectDiscreteGenerators(DiscrP, DiscrZ);
    for (int iG = 0; iG < DiscrGen.size(); iG++) {
      cctbx_assert(m_VM.size() < 3);
      TrVec G = DiscrGen[iG].cancel();
      ssVM VM;
      VM.V = G.vec();
      VM.M = G.BF();
      m_VM.push_back(VM);
    }
    std::sort(m_VM.begin(), m_VM.end(), Cmp_ssVM());
  }

  bool StructureSeminvariant::is_ss(const miller::Index& H) const
  {
    for(std::size_t i=0;i<m_VM.size();i++) {
      int u = af::sum(m_VM[i].V * H);
      if (m_VM[i].M) {
        if (u % m_VM[i].M) return false;
      }
      else {
        if (u)             return false;
      }
    }
    return true;
  }

  af::small<int, 3>
  StructureSeminvariant::apply_mod(const miller::Index& H) const
  {
    af::small<int, 3> result;
    for(std::size_t i=0;i<m_VM.size();i++) {
      result.push_back(af::sum(m_VM[i].V * H));
      if (m_VM[i].M) result[i] %= m_VM[i].M;
    }
    return result;
  }

}} // namespace cctbx::sgtbx
