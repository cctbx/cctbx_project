// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 02: start port of sglite/sgss.c (R.W. Grosse-Kunstleve)
 */

#include <algorithm>
#include <cctbx/sgtbx/select_generators.h>
#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/reference.h>

namespace sgtbx {

  namespace seminvariant_cpp { // anonymous namespace causes
                               // warnings with some compilers

    inline void copy(const int *source, int* target, std::size_t n) {
      for(std::size_t i=0;i<n;i++) target[i] = source[i];
    }

    boost::array<int, 3 * 3*3>
    ConstructGenRmI(const detail::AnyGenerators& Gen, bool Primitive)
    {
      boost::array<int, 3 * 3*3> result;
      for(std::size_t i=0;i<Gen.nGen;i++) {
        const RotMx R = Gen.ZGen[i].Rpart();
        if (!Primitive) {
          copy(R.minusUnit().elems, &result.elems[i * 3*3], 3*3);
        }
        else {
          copy(Gen.Z2POp(R).minusUnit().elems, &result.elems[i * 3*3], 3*3);
        }
      }
      if (Gen.ZInvT.isValid()) {
        copy(RotMx(1, -2).elems, &result.elems[Gen.nGen * 3*3], 3*3);
      }
      return result;
    }

    class CmpOLen2 {
      public:
        CmpOLen2() : CmpT(3) {}
        bool operator()(const Vec3& a, const Vec3& b) const
        {
          int i;
          int OLen2a = 0; for(i=0;i<3;i++) OLen2a += a[i] * a[i];
          int OLen2b = 0; for(i=0;i<3;i++) OLen2b += b[i] * b[i];
          if (OLen2a < OLen2b) return true;
          if (OLen2a > OLen2b) return false;
          return CmpT(a.elems, b.elems);
        }
      private:
        const CmpiVect CmpT;
    };

    cctbx::fixcap_vector<ssVM, 3>
    GetContNullSpace(const detail::AnyGenerators& Gen)
    {
      cctbx::fixcap_vector<ssVM, 3> result;
      boost::array<int, 3 * 3 * 3> GenRmI = ConstructGenRmI(Gen, false);
      int RankGenRmI = iRowEchelonFormT(GenRmI.elems, Gen.nAll() * 3, 3, 0, 0);
      cctbx_assert(RankGenRmI >= 0 && RankGenRmI <= 3);
      int IxIndep[3];
      int nIndep = iRESetIxIndep(GenRmI.elems, RankGenRmI, 3, IxIndep, 3);
      cctbx_assert(nIndep >= 0);
      if (nIndep != 2) {
        for (int iIndep = 0; iIndep < nIndep; iIndep++) {
          ssVM VM;
          if (Gen.nAll() == 0) VM.V.assign(0);
          VM.V[IxIndep[iIndep]] = 1;
          cctbx_assert(iREBacksubst(GenRmI.elems, 0, RankGenRmI, 3,
                                    VM.V.elems, 0) > 0);
          VM.M = 0;
          result.push_back(VM);
        }
      }
      else {
        boost::array<Vec3, 4> Sol;
        SolveHomRE1(GenRmI.elems, IxIndep, Sol.elems);
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

    typedef cctbx::fixcap_vector<TrVec, 8> DiscrList;

    void UpdateBestZ(const DiscrList& OrigZf,
                     const TrVec& Shift,
                     DiscrList& BestZf,
                     DiscrList& BestZc)
    {
      for (int iDL = 1; iDL < OrigZf.size(); iDL++) {
        TrVec Zf = (OrigZf[iDL] + Shift).modPositive();
        TrVec Zc = Zf.cancel();
        for(int i=0;i<3;i++) {
          if (Zf[i]) {
            if (   CmpOLen2()(Zc, BestZc[iDL])
                || (   static_cast<Vec3>(Zc) == static_cast<Vec3>(BestZc[iDL])
                    && Zc.BF() < BestZc[iDL].BF())) {
              BestZf[iDL] = Zf;
              BestZc[iDL] = Zc;
            }
            break;
          }
        }
      }
    }

    void BestVectors(const SgOps& sgo,
                     const cctbx::fixcap_vector<ssVM, 3>& ContinuousVM,
                     DiscrList& DiscrZ)
    {
      if (sgo.nLTr() == 1 && ContinuousVM.size() == 0) return;
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
      for (iLTr = 1; iLTr < sgo.nLTr(); iLTr++) {
        int BF = sgo.LTr(iLTr).BF();
        for(int i=0;i<3;i++) {
          int g = gcd(sgo.LTr(iLTr)[i], BF);
          LTBF = lcm(LTBF, BF / g);
        }
      }
      int fGrd = 1;
      int iVM;
      for (iVM = 0; iVM < ContinuousVM.size(); iVM++) {
        for(int i=0;i<3;i++) {
          if (ContinuousVM[iVM].V[i]) fGrd = lcm(fGrd, ContinuousVM[iVM].V[i]);
        }
      }
      LTBF *= fGrd;
      cctbx_assert(LTBF > 0);
      if (LTBF > 6) LTBF = lcm(LTBF,  6);
      else          LTBF = lcm(LTBF, 12);
      DiscrList OrigZf;
      DiscrList BestZf;
      DiscrList BestZc;
      for (iDL = 0; iDL < DiscrZ.size(); iDL++) {
        OrigZf.push_back(DiscrZ[iDL].newBaseFactor(LTBF));
        BestZf.push_back(OrigZf[iDL]);
        BestZc.push_back(BestZf[iDL].cancel());
      }
      cctbx_assert(ContinuousVM.size() < 3);
      cctbx::fixcap_vector<int, 2> loop_min, loop_max;
      for (iVM = 0; iVM < ContinuousVM.size(); iVM++) {
        loop_min.push_back(0);
        loop_max.push_back(LTBF - 1);
      }
      TrVec LTr[3];
      for (iLTr = 0; iLTr < sgo.nLTr(); iLTr++) {
        LTr[0] = sgo.LTr(iLTr).newBaseFactor(LTBF);
        NestedLoop<cctbx::fixcap_vector<int, 2> > loop(loop_min, loop_max);
        do {
          for (iVM = 0; iVM < ContinuousVM.size(); iVM++) {
            const cctbx::fixcap_vector<int, 2>& f = loop();
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

    cctbx::fixcap_vector<TrVec, 3>
    SelectDiscreteGenerators(const DiscrList& DiscrP,
                             const DiscrList& DiscrZ)
    {
      if (DiscrP.size() == 1) return cctbx::fixcap_vector<TrVec, 3>();
      for (int nIx = 1; nIx <= DiscrP.size() - 1 && nIx <= 3; nIx++) {
        int Ix[3], iIx;
        for (iIx = 0; iIx < nIx; iIx++) Ix[iIx] = iIx;
        do {
          cctbx::fixcap_vector<TrVec, 3> result;
          TrOps DiscrGrp(DiscrP[0].BF());
          for (iIx = 0; iIx < nIx; iIx++) {
            DiscrGrp.expand(DiscrP[Ix[iIx] + 1]);
            result.push_back(DiscrZ[Ix[iIx] + 1]);
          }
          if (DiscrGrp.nVects() == DiscrP.size()) return result;
          cctbx_assert(DiscrGrp.nVects() < DiscrP.size());
        }
        while (NextOf_n_from_m(DiscrP.size() - 1, nIx, Ix) != 0);
      }
      throw cctbx_internal_error();
    }

    class CmpDiscrZ {
      public:
        CmpDiscrZ() : CmpT(3) {}
        inline bool operator()(const TrVec& a,
                               const TrVec& b) const {
          return CmpT(a.elems, b.elems);
        }
      private:
        const CmpiVect CmpT;
    };

    class Cmp_ssVM {
      public:
        Cmp_ssVM() : CmpT(3) {}
        inline bool operator()(const ssVM& a, const ssVM& b) const {
          return CmpT(a.V.elems, b.V.elems);
        }
      private:
        const CmpiVect CmpT;
    };

  } // namespace seminvariant_cpp

  using namespace seminvariant_cpp;

  StructureSeminvariant::StructureSeminvariant(const SgOps& sgo)
  {
    detail::AnyGenerators Gen(sgo);
    m_VM = GetContNullSpace(Gen);
    if (m_VM.size() == 3) return; // space group P1
    boost::array<int, 3 * 3 * 3> SNF = ConstructGenRmI(Gen, true);
    int Q[3 * 3];
    int nd = SmithNormalForm(SNF.elems, Gen.nAll() * 3, 3, 0, Q);
    cctbx_assert(nd >=0 && nd <= 3);
    int id;
    int DTBF = 1;
    for (id = 0; id < nd; id++) DTBF = lcm(DTBF, SNF[(nd + 1) * id]);
    TrOps DiscrGrpP(DTBF);
    for (id = 0; id < nd; id++) {
      int d = SNF[(nd + 1) * id];
      for (int f = 1; f < d; f++) {
        Vec3 xp;
        xp.assign(0);
        xp[id] = f * DTBF / d;
        TrVec x(DTBF);
        MatrixLite::multiply<int>(Q, xp.elems, 3, 3, 1, x.elems);
        DiscrGrpP.expand(x);
      }
    }
    int iDL;
    DiscrList DiscrZ;
    for (iDL = 0; iDL < DiscrGrpP.nVects(); iDL++) {
      DiscrZ.push_back(
        (Gen.Z2POp.InvM().Rpart() * DiscrGrpP[iDL]).modPositive());
    }
    BestVectors(sgo, m_VM, DiscrZ);
    std::sort(DiscrZ.begin(), DiscrZ.end(), CmpDiscrZ());
    DiscrList DiscrP;
    for (iDL = 0; iDL < DiscrZ.size(); iDL++) {
      DiscrP.push_back(
        (Gen.Z2POp.M().Rpart() * DiscrZ[iDL]).newBaseFactor(DTBF));
    }
    cctbx::fixcap_vector<TrVec, 3>
    DiscrGen = SelectDiscreteGenerators(DiscrP, DiscrZ);
    for (int iG = 0; iG < DiscrGen.size(); iG++) {
      cctbx_assert(m_VM.size() < 3);
      TrVec G = DiscrGen[iG].cancel();
      ssVM VM;
      VM.V = static_cast<Vec3>(G);
      VM.M = G.BF();
      m_VM.push_back(VM);
    }
    std::sort(m_VM.begin(), m_VM.end(), Cmp_ssVM());
  }

  bool StructureSeminvariant::is_ss(const Miller::Index& H) const
  {
    for(std::size_t i=0;i<m_VM.size();i++) {
      int u = m_VM[i].V * H;
      if (m_VM[i].M) {
        if (u % m_VM[i].M) return false;
      }
      else {
        if (u)             return false;
      }
    }
    return true;
  }

  cctbx::fixcap_vector<int, 3>
  StructureSeminvariant::apply_mod(const Miller::Index& H) const
  {
    cctbx::fixcap_vector<int, 3> result;
    for(std::size_t i=0;i<m_VM.size();i++) {
      result.push_back(m_VM[i].V * H);
      if (m_VM[i].M) result[i] %= m_VM[i].M;
    }
    return result;
  }

} // namespace sgtbx
