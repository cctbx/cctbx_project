// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 02: start port of sglite/sgss.c (R.W. Grosse-Kunstleve)
 */

#include <algorithm>
#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/reference.h>

namespace sgtbx {

  namespace detail {

    struct DiscrList { // XXX split this up
      TrVec P; // primitive cell
      TrVec Z; // centred cell
    };

    struct AnyGenerators { // see also: StdGenerators
      AnyGenerators(const SgOps& sgo);
      inline int nAll() const {
        if (ZInvT.isValid()) return nGen + 1;
        return nGen;
      }
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

    void copy(const int *source, int* target, std::size_t n) {
      for(std::size_t i=0;i<n;i++) target[i] = source[i];
    }

  } // namespace detail

  boost::array<int, 3 * 3 * 3>
  ConstructGenRmI(const detail::AnyGenerators& Gen, bool Primitive)
  {
    boost::array<int, 3 * 3 * 3> result;
    for(std::size_t i=0;i<Gen.nGen;i++) {
      const RotMx R = Gen.ZGen[i].Rpart();
      if (!Primitive) {
        detail::copy(R.minusUnit().elems,
                     &result.elems[i * 3 * 3], 3 * 3);
      }
      else {
        detail::copy(Gen.Z2POp(R).minusUnit().elems,
                     &result.elems[i * 3 * 3], 3 * 3);
      }
    }
    if (Gen.ZInvT.isValid()) {
      detail::copy(RotMx(1, -2).elems,
                   &result.elems[Gen.nGen * 3 * 3], 3 * 3);
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

  void
  StructureSeminvariant::GetContNullSpace(const detail::AnyGenerators& Gen)
  {
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
        m_VM.push_back(VM);
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
        m_VM.push_back(VM);
      }
    }
  }

  void UpdateBestZ(int nDiscrLst,
                   const TrVec OrigZf[detail::mDiscrList],
                   const TrVec& Shift,
                   TrVec BestZf[detail::mDiscrList],
                   TrVec BestZc[detail::mDiscrList])
  {
    for (int iDL = 1; iDL < nDiscrLst; iDL++) {
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

  void StructureSeminvariant::BestVectors(
    const SgOps& sgo,
    cctbx::static_vector<detail::DiscrList,
                         detail::mDiscrList>& DiscrLst) const
  {
    if (sgo.nLTr() == 1 && m_VM.size() == 0) return;
    int LTBF = 1;
    int iDL;
    for (iDL = 1; iDL < DiscrLst.size(); iDL++) {
      int BF = DiscrLst[iDL].Z.BF();
      for(int i=0;i<3;i++) {
        int g = gcd(DiscrLst[iDL].Z[i], BF);
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
    for (iVM = 0; iVM < m_VM.size(); iVM++) {
      for(int i=0;i<3;i++) {
        if (m_VM[iVM].V[i]) fGrd = lcm(fGrd, m_VM[iVM].V[i]);
      }
    }
    LTBF *= fGrd;
    cctbx_assert(LTBF > 0);
    if (LTBF > 6) LTBF = lcm(LTBF,  6);
    else          LTBF = lcm(LTBF, 12);
    TrVec OrigZf[detail::mDiscrList];
    TrVec BestZf[detail::mDiscrList];
    TrVec BestZc[detail::mDiscrList];
    for (iDL = 1; iDL < DiscrLst.size(); iDL++) {
      OrigZf[iDL] = DiscrLst[iDL].Z.newBaseFactor(LTBF);
      BestZf[iDL] = OrigZf[iDL];
      BestZc[iDL] = BestZf[iDL].cancel();
    }
    cctbx_assert(m_VM.size() < 3);
    cctbx::static_vector<int, 2> loop_min, loop_max;
    for (iVM = 0; iVM < m_VM.size(); iVM++) {
      loop_min.push_back(0);
      loop_max.push_back(LTBF - 1);
    }
    TrVec LTr[3];
    for (iLTr = 0; iLTr < sgo.nLTr(); iLTr++) {
      LTr[0] = sgo.LTr(iLTr).newBaseFactor(LTBF);
      NestedLoop<cctbx::static_vector<int, 2> > loop(loop_min, loop_max);
      do {
        for (iVM = 0; iVM < m_VM.size(); iVM++) {
          const cctbx::static_vector<int, 2>& f = loop();
          LTr[iVM + 1] = LTr[iVM] + f[iVM] * TrVec(m_VM[iVM].V, LTBF);
        }
        UpdateBestZ(DiscrLst.size(), OrigZf, LTr[m_VM.size()], BestZf, BestZc);
        loop.incr();
      }
      while (loop.over() == 0);
    }
    for (iDL = 1; iDL < DiscrLst.size(); iDL++) {
      DiscrLst[iDL].Z = BestZf[iDL].newBaseFactor(DiscrLst[iDL].Z.BF());
    }
  }

  cctbx::static_vector<TrVec, 3>
  SelectDiscreteGenerators(cctbx::static_vector<detail::DiscrList,
                                                detail::mDiscrList> DiscrLst)
  {
    if (DiscrLst.size() == 1) return cctbx::static_vector<TrVec, 3>();
    for (int nIx = 1; nIx <= DiscrLst.size() - 1 && nIx <= 3; nIx++) {
      int Ix[3], iIx;
      for (iIx = 0; iIx < nIx; iIx++) Ix[iIx] = iIx;
      do {
        cctbx::static_vector<TrVec, 3> result;
        TrOps DiscrGrp(DiscrLst[0].P.BF());
        for (iIx = 0; iIx < nIx; iIx++) {
          DiscrGrp.expand(DiscrLst[Ix[iIx] + 1].P);
          result.push_back(DiscrLst[Ix[iIx] + 1].Z);
        }
        if (DiscrGrp.nVects() == DiscrLst.size()) return result;
        cctbx_assert(DiscrGrp.nVects() < DiscrLst.size());
      }
      while (NextOf_n_from_m(DiscrLst.size() - 1, nIx, Ix) != 0);
    }
    throw cctbx_internal_error();
  }

  class CmpDiscr {
    public:
      CmpDiscr() : CmpT(3) {}
      inline bool operator()(const detail::DiscrList& a,
                             const detail::DiscrList& b) const {
        return CmpT(a.Z.elems, b.Z.elems);
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

  StructureSeminvariant::StructureSeminvariant(const SgOps& sgo)
  {
    detail::AnyGenerators Gen(sgo);
    GetContNullSpace(Gen);
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
    cctbx_assert(DiscrGrpP.nVects() <= detail::mDiscrList);
    cctbx::static_vector<detail::DiscrList, detail::mDiscrList> DiscrLst;
    for (int iDL = 0; iDL < DiscrGrpP.nVects(); iDL++) {
      detail::DiscrList v;
      v.P = DiscrGrpP[iDL];
      v.Z = (Gen.Z2POp.InvM().Rpart() * v.P).modPositive();
      DiscrLst.push_back(v);
    }
    cctbx_assert(DiscrLst.size() == DiscrGrpP.nVects());
    BestVectors(sgo, DiscrLst);
    std::sort(DiscrLst.begin(), DiscrLst.end(), CmpDiscr());
    cctbx::static_vector<TrVec, 3>
    DiscrGen = SelectDiscreteGenerators(DiscrLst);
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

  Vec3 StructureSeminvariant::get_uvw(const Miller::Index& H) const
  {
    Vec3 result;
    for(std::size_t i=0;i<m_VM.size();i++) {
      result[i] = m_VM[i].V * H;
      if (m_VM[i].M) result[i] %= m_VM[i].M;
    }
    return result;
  }

} // namespace sgtbx
