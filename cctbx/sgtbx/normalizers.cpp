// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Created: 10-May-2001 (R.W. Grosse-Kunstleve)
 */

#include <boost/rational.hpp>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/reference.h>

namespace cctbx { namespace sgtbx {
  namespace ReferenceSettings {

    std::vector<RTMx>
    GetNormAddlG(int SgNumber, bool affine, bool UseK2L, bool UseL2N)
    {
      using sgtbx::tables::ReferenceSettings::NormalizerAddlGenerators;

      cctbx_assert(0 < SgNumber && SgNumber <= 230);
      std::vector<RTMx> result;
      for (int iType = 0; iType < 2; iType++) {
        const char* HallMxSymbol = 0;
        if      (iType == 0 && UseK2L)
          HallMxSymbol = NormalizerAddlGenerators[SgNumber].K2L;
        else if (iType == 1 && UseL2N && (SgNumber >= 75 || affine))
          HallMxSymbol = NormalizerAddlGenerators[SgNumber].L2N;
        if (HallMxSymbol) {
          SpaceGroup SgOpsAddlG(true); // NoExpand: no group multiplication
          parse_string ps(HallMxSymbol);
          int nAddedMx = SgOpsAddlG.ParseHallSymbol(ps, true, true);
          cctbx_assert(nAddedMx > 0);
          cctbx_assert(SgOpsAddlG.nLTr() == 1);
          if (SgOpsAddlG.isCentric()) {
            result.push_back(SgOpsAddlG(0, 1, 0));
          }
          for (int i = 1; i < SgOpsAddlG.nSMx(); i++) {
            result.push_back(SgOpsAddlG[i]);
          }
        }
      }
      return result;
    }

    void GetMonoAffNormTrialRanges(const RotMx& CBMxR, int& r00, int& r22)
    {
      /* International Tables Volume A, chapter 15, tables 15.3.3 & 15.3.4.

       M.C = n00*c00 + n02*c20,  n00*c01 + n02*c21,  n00*c02 + n02*c22,
             c10,                c11,                c12,
             n20*c00 + n22*c20,  n20*c01 + n22*c21,  n20*c02 + n22*c22

       Determine trial range for n00 and n20:
         max(lcm(c00, c20) / c00,
             lcm(c01, c21) / c01,
             lcm(c02, c22) / c02)
       Determine trial range for n02 and n22:
         max(lcm(c00, c20) / c20,
             lcm(c01, c21) / c21,
             lcm(c02, c22) / c22)
       */

      r00 = 1;
      r22 = 1;
      for(int i=0;i<3;i++) {
        int l = boost::lcm(CBMxR[i], CBMxR[6 + i]);
        if (CBMxR[i]) {
          int n = std::abs(l / CBMxR[i]);
          if (r00 < n) r00 = n;
        }
        if (CBMxR[i + 6]) {
          int n = std::abs(l / CBMxR[6 + i]);
          if (r22 < n) r22 = n;
        }
      }
      r00++;
      r22++;
    }

    bool CheckMonoAffNormRestrictions(int SgNumber, const RotMx& M)
    {
      // International Tables Volume A, chapter 15, tables 15.3.3 & 15.3.4.
      switch (SgNumber) {
        case  3:
        case  4:
        case  6:
        case 10:
        case 11: /* M2 */
          break;

        case  5:
        case  8:
        case 12: /* M4 */
        case  9:
        case 15: /* M6 or M12 */
          if (M[0] % (2 * M.BF()) == 0) return false;
          if (M[6] % (2 * M.BF()) != 0) return false;
          if (M[8] % (2 * M.BF()) == 0) return false;
          break;

        case  7:
        case 13:
        case 14: /* M5 */
          if (M[0] % (2 * M.BF()) == 0) return false;
          if (M[2] % (2 * M.BF()) != 0) return false;
          if (M[8] % (2 * M.BF()) == 0) return false;
          break;

        default:
          throw cctbx_internal_error();
      }
      return true;
    }

  } // namespace ReferenceSettings

  std::vector<RTMx>
  SpaceGroupInfo::getAddlGeneratorsOfEuclideanNormalizer(bool getK2L,
                                                         bool getL2N) const
  {
    std::vector<RTMx> result
      = ReferenceSettings::GetNormAddlG(SgNumber(), false, getK2L, getL2N);
    ChOfBasisOp C = CBOp().swap();
    for(std::size_t i=0;i<result.size();i++) {
      result[i] = C(result[i]);
    }
    return result;
  }

}} // namespace cctbx::sgtbx
