// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 10: created from fragments in type.cpp, seminvariant.cpp (rwgk)
 */

#include <cctbx/sgtbx/select_generators.h>

namespace cctbx { namespace sgtbx {
  namespace detail {

    AnyGenerators::AnyGenerators(const SpaceGroup& SgOps) : nGen(0)
    {
      using namespace tables::CrystalSystem;

      Z2POp = SgOps.getZ2POp();

      ZInvT = SgOps.InvT(true);
      PInvT = TrVec(0);
      int i;
      for(i=0;i<2;i++) ZGen[i] = RTMx(0, 0);
      for(i=0;i<2;i++) PGen[i] = RTMx(0, 0);

      int PrincipalProperOrder = 0;

      tables::MatrixGroup::Code PG_MGC = SgOps.getPointGroupType();
      switch (PG_MGC.CrystalSystem())
      {
        case Triclinic:
          break;

        case Monoclinic:
          ZGen[0] = SgOps[1];
          nGen = 1;
          break;

        case Orthorhombic:
          ZGen[0] = SgOps[1];
          ZGen[1] = SgOps[2];
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
            for(i = 1; i < SgOps.nSMx(); i++) {
              PrincipalRI = SgOps[i].Rpart().getInfo();
              if (std::abs(PrincipalRI.Rtype()) == PrincipalProperOrder) {
                if (PrincipalRI.SenseOfRotation() > 0) {
                  ZGen[0] = SgOps[i];
                  nGen++;
                  break;
                }
              }
            }
            cctbx_assert(nGen == 1);
            int iPrincipal = i;
            for(i = 1; i < SgOps.nSMx(); i++) {
              if (i == iPrincipal) continue;
              RotMxInfo RI = SgOps[i].Rpart().getInfo();
              if (std::abs(RI.Rtype()) == 2) {
                if (PrincipalRI.EV() != RI.EV()) {
                  ZGen[1] = SgOps[i];
                  nGen++;
                  break;
                }
              }
            }
          }
          break;

        case Cubic:
          for(i = 1; i < SgOps.nSMx(); i++) {
            RotMxInfo RI = SgOps[i].Rpart().getInfo();
            if      (std::abs(RI.Rtype()) == 3) {
              if (RI.SenseOfRotation() > 0) {
                if (!ZGen[0].isValid()) {
                  ZGen[0] = SgOps[i];
                  nGen++;
                  if (nGen == 2) break;
                }
              }
            }
            else if (std::abs(RI.Rtype()) == SgOps.nSMx() / 6) {
              if (RI.SenseOfRotation() >= 0) {
                if (!ZGen[1].isValid()) {
                  ZGen[1] = SgOps[i];
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

    StdGenerators::StdGenerators(const SpaceGroup& WorkSgOps,
                                 const tables::MatrixGroup::Code& PG_MGC)
      : AnyGenerators()
    {
      using namespace tables::CrystalSystem;

      const af::int3 EV_001( 0, 0, 1);
      const af::int3 EV_100( 1, 0, 0);
      const af::int3 EV_110( 1, 1, 0);
      const af::int3 EV_m10(-1, 1, 0);
      const af::int3 EV_111( 1, 1, 1);

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
      for (i = 0; i < nGen; i++) ZGen[i].modPositiveInPlace();
    }

  } // namespace detail
}} // namespace cctbx::sgtbx
