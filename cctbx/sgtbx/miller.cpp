// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 31: Redesign: AsymIndex (rwgk)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_ops.h>
#include <cctbx/sgtbx/miller_ref_asu.h>

namespace cctbx { namespace sgtbx {

  PhaseInfo::PhaseInfo(
    SpaceGroup const& sgops,
    miller::Index const& h,
    bool no_test_sys_absent)
    : m_HT(-1), m_TBF(sgops.TBF()), m_SysAbsWasTested(!no_test_sys_absent)
  {
    // systematically absent reflection: if HR == H and HT != 0 mod 1
    // restricted phase: if HR == -H: phi(H) = pi*HT + n*pi
    if (no_test_sys_absent) {
      // Fast determination of phase restriction without considering
      // conditions for systematically absent reflections.
      if (sgops.isCentric()) {
        m_HT = HT_mod_1(h, sgops.InvT());
      }
      else {
        for(int i=0;i<sgops.nSMx();i++) {
          if (h * sgops[i].Rpart() == -h) {
            m_HT = HT_mod_1(h, sgops[i].Tpart());
            break;
          }
        }
      }
      return;
    }
    // Simulatenous determination of phase restriction and evaluation
    // conditions for systematically absent reflections.
    for(int iSMx=0;iSMx<sgops.nSMx();iSMx++) {
      const RotMx& R = sgops[iSMx].Rpart();
      const TrVec& T = sgops[iSMx].Tpart();
      TrVec TS(0);
      TrVec TR(0);
      miller::Index HR = h * R;
      if      (h == HR) {
        TS = T;
        if (sgops.isCentric()) TR = sgops.InvT() - T;
      }
      else if (h == -HR) {
        TR = T;
        if (sgops.isCentric()) TS = sgops.InvT() - T;
      }
      if (TS.isValid()) {
        for(int iLTr=0;iLTr<sgops.nLTr();iLTr++) {
          if ((h * (TS + sgops.LTr(iLTr))) % TS.BF() != 0) {
            m_HT = -2;
            return;
          }
        }
      }
      if (TR.isValid()) {
        for(int iLTr=0;iLTr<sgops.nLTr();iLTr++) {
          int HT = HT_mod_1(h, TR + sgops.LTr(iLTr));
          if      (m_HT < 0) m_HT = HT;
          else if (m_HT != HT) {
            m_HT = -2;
            return;
          }
        }
      }
    }
  }

  bool
  PhaseInfo::isValidPhase(double phi, bool deg, double tolerance) const
  {
    if (!isCentric()) return true;
    double period = ht_period(deg);
    double delta = std::fmod(phi - HT_angle(deg), period);
    if (delta >  tolerance) delta -= period;
    if (delta < -tolerance) delta += period;
    if (delta <= tolerance) return true;
    return false;
  }

  double
  PhaseInfo::nearest_valid_phase(double phi, bool deg) const
  {
    if (!isCentric()) return phi;
    double period = ht_period(deg);
    double phi_restr = HT_angle(deg);
    double delta0 = phi - phi_restr;
    double delta1 = delta0 - period;
    if (  math::abs(std::fmod(delta1, period))
        < math::abs(std::fmod(delta0, period))) return phi_restr + period;
    return phi_restr;
  }

  af::shared<bool>
  SpaceGroup::isSysAbsent(af::shared<miller::Index> H) const
  {
    af::shared<bool> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(isSysAbsent(H[i]));
    }
    return result;
  }

  bool SpaceGroup::isCentric(const miller::Index& H) const
  {
    if (isCentric()) return true;
    for(std::size_t i=0;i<m_nSMx;i++) {
      if (H * m_SMx[i].Rpart() == -H) return true;
    }
    return false;
  }

  af::shared<bool>
  SpaceGroup::isCentric(const af::shared<miller::Index>& H) const
  {
    af::shared<bool> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(isCentric(H[i]));
    }
    return result;
  }

  int SpaceGroup::epsilon(const miller::Index& H) const
  {
    int result = 0;
    for(std::size_t i=0;i<m_nSMx;i++) {
      miller::Index HR = H * m_SMx[i].Rpart();
      if (HR == H || (isCentric() && HR == -H))
        result++;
    }
    cctbx_assert(result != 0 && m_nSMx % result == 0);
    return result;
  }

  af::shared<int>
  SpaceGroup::epsilon(const af::shared<miller::Index>& H) const
  {
    af::shared<int> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(epsilon(H[i]));
    }
    return result;
  }

  int SpaceGroup::multiplicity(const miller::Index& H, bool FriedelFlag) const
  {
    if (H.is000()) return 1;
    int Centro = (isCentric() || FriedelFlag);
    int M = 0;
    int R = 0;
    for(std::size_t i=0;i<m_nSMx;i++) {
      miller::Index HR = H * m_SMx[i].Rpart();
      if      (HR == H) M++;
      else if (HR == -H) R++;
    }
    cctbx_assert(M != 0 && m_nSMx % M == 0 && (R == 0 || R == M));
    M = m_nSMx / M;
    if (Centro && R == 0) M *= 2;
    return M;
  }

  af::shared<int>
  SpaceGroup::multiplicity(const af::shared<miller::Index>& H,
                           bool FriedelFlag) const
  {
    af::shared<int> result;
    result.reserve(H.size());
    for(std::size_t i=0;i<H.size();i++) {
      result.push_back(multiplicity(H[i], FriedelFlag));
    }
    return result;
  }

}} // namespace cctbx::sgtbx
