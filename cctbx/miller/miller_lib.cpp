// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jul 2002: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller/build.h>
#include <cctbx/miller/join.h>
#include <cctbx/miller/bins.h>

namespace cctbx { namespace miller {

  SymEquivIndices::SymEquivIndices(
    sgtbx::SpaceGroup const& sgops,
    Index const& h_in)
  : m_TBF(sgops.TBF()),
    m_OrderP(sgops.OrderP()),
    m_HT_Restriction(-1)
  {
    using namespace sgtbx;
    int iInv;
    for(iInv=0;iInv<sgops.fInv();iInv++) {
      int iSMx;
      for(iSMx=0;iSMx<sgops.nSMx();iSMx++) {
        RTMx M = sgops(0, iInv, iSMx);
        Index HR = h_in * M.Rpart();
        bool found = false;
        for(int i=0;i<N();i++) {
          if (m_List[i].HR() == HR) {
            found = true;
            break;
          }
        }
        if (!found) {
          add(SymEquivIndex(HR, HT_mod_1(h_in, M.Tpart()), m_TBF, false));
        }
      }
    }
    cctbx_assert((sgops.nSMx() * sgops.fInv()) % N() == 0);
    cctbx_assert(!isCentric() || N() % 2 == 0);
  }

  void SymEquivIndices::add(SymEquivIndex const& SEI)
  {
    m_List.push_back(SEI);
    if (m_List.size() > 1) {
      if (SEI.HR() == -m_List[0].HR()) {
        cctbx_assert(m_HT_Restriction < 0 || m_HT_Restriction == SEI.HT());
        m_HT_Restriction = SEI.HT();
      }
    }
  }

  SymEquivIndex
  SymEquivIndices::operator()(int iMate, int iList) const
  {
    if (   iMate < 0 || iMate >= fMates(true)
        || iList < 0 || iList >= N()) {
      throw error_index();
    }
    return m_List[iList].Mate(iMate);
  }

  SymEquivIndices::iIL_decomposition
  SymEquivIndices::decompose_iIL(int iIL) const
  {
    // iIL = iMate * N + iList
    if (iIL < 0 || iIL >= M(true)) {
      throw error_index();
    }
    return iIL_decomposition(iIL / N(), iIL % N());
  }

  SymEquivIndex
  SymEquivIndices::operator()(int iIL) const
  {
    iIL_decomposition d = decompose_iIL(iIL);
    return operator()(d.iMate, d.iList);
  }

  af::shared<SymEquivIndex>
  SymEquivIndices::p1_listing(bool friedel_flag) const
  {
    af::shared<SymEquivIndex> result;
    if (!friedel_flag) {
      result.reserve(N());
      for(std::size_t i=0;i<N();i++) result.push_back(m_List[i]);
    }
    else {
      if (isCentric()) result.reserve(N() / 2);
      else             result.reserve(N());
      for(std::size_t i=0;i<M(true);i++) {
        SymEquivIndex h_eq = operator()(i);
        if (sgtbx::isInReferenceReciprocalSpaceASU_1b(h_eq.H())) {
          result.push_back(h_eq);
        }
      }
    }
    cctbx_assert(result.size() == result.capacity());
    return result;
  }

  void IndexGenerator::InitializeLoop(Index const& ReferenceHmax)
  {
    af::int3 CutP = m_ASU.ReferenceASU()->getCutParameters();
    Index ReferenceHbegin;
    Index ReferenceHend;
    for(std::size_t i=0;i<3;i++) {
      ReferenceHbegin[i] = ReferenceHmax[i] * CutP[i];
      ReferenceHend[i] = ReferenceHmax[i] + 1;
    }
    m_loop = af::nested_loop<Index>(ReferenceHbegin, ReferenceHend);
    m_next_is_minus_previous = false;
  }

  IndexGenerator::IndexGenerator(uctbx::UnitCell const& uc,
                                 sgtbx::SpaceGroupInfo const& SgInfo,
                                 bool FriedelFlag,
                                 double Resolution_d_min)
    : m_UnitCell(uc),
      m_SgNumber(SgInfo.SgNumber()),
      m_SgOps(SgInfo.SgOps()),
      m_FriedelFlag(FriedelFlag),
      m_ASU(SgInfo)
  {
    if (Resolution_d_min <= 0.) {
      throw error("Resolution limit must be greater than zero.");
    }
    m_Qhigh = 1. / (Resolution_d_min * Resolution_d_min);
    uctbx::UnitCell
    ReferenceUnitCell = m_UnitCell.ChangeBasis(SgInfo.CBOp().InvM().Rpart());
    InitializeLoop(ReferenceUnitCell.MaxMillerIndices(Resolution_d_min));
  }

  IndexGenerator::IndexGenerator(sgtbx::SpaceGroupInfo const& SgInfo,
                                 bool FriedelFlag,
                                 Index const& MaxIndex)
    : m_UnitCell(),
      m_SgNumber(SgInfo.SgNumber()),
      m_SgOps(SgInfo.SgOps()),
      m_FriedelFlag(FriedelFlag),
      m_ASU(SgInfo),
      m_Qhigh(-1.)
  {
    InitializeLoop(Index(af::abs(MaxIndex)));
  }

  bool IndexGenerator::set_phase_info(Index const& h)
  {
    m_phase_info = sgtbx::PhaseInfo(m_SgOps, h, false);
    return m_phase_info.isSysAbsent();
  }

  Index IndexGenerator::next_under_friedel_symmetry()
  {
    const int RBF = m_ASU.CBOp().M().RBF();
    for (; m_loop.over() == 0;) {
      Index ReferenceH = m_loop();
      m_loop.incr();
      if (m_ASU.ReferenceASU()->isInASU(ReferenceH)) {
        if (m_ASU.isReferenceASU()) {
          if (m_Qhigh < 0.) {
            if (!ReferenceH.is000() && !set_phase_info(ReferenceH)) {
              return ReferenceH;
            }
          }
          else {
            double Q = m_UnitCell.Q(ReferenceH);
            if (Q != 0 && Q <= m_Qhigh && !set_phase_info(ReferenceH)) {
              return ReferenceH;
            }
          }
        }
        else {
          sgtbx::TrVec HR(ReferenceH * m_ASU.CBOp().M().Rpart(), RBF);
          HR = HR.cancel();
          if (HR.BF() == 1) {
            Index H(HR.vec());
            if (m_Qhigh < 0.) {
              if (!H.is000() && !set_phase_info(H)) {
                return H;
              }
            }
            else {
              double Q = m_UnitCell.Q(H);
              if (Q != 0 && Q <= m_Qhigh && !set_phase_info(H)) {
                return H;
              }
            }
          }
        }
      }
    }
    return Index(0, 0, 0);
  }

  Index IndexGenerator::next()
  {
    if (m_FriedelFlag) return next_under_friedel_symmetry();
    if (m_next_is_minus_previous) {
      m_next_is_minus_previous = false;
      return -m_previous;
    }
    m_previous = next_under_friedel_symmetry();
    if (m_previous.is000()) return m_previous;
    m_next_is_minus_previous = !m_phase_info.isCentric();
    if (m_next_is_minus_previous && 143 <= m_SgNumber && m_SgNumber <= 167) {
      // For trigonal space groups it has to be checked if a symmetrically
      // equivalent index of the Friedel opposite is in the ASU.
      cctbx_assert(!m_SgOps.isCentric());
      Index minus_h = -m_previous;
      for(int i=0;i<m_SgOps.nSMx();i++) {
        Index minus_h_eq = minus_h * m_SgOps[i].Rpart();
        if (m_ASU.isInASU(minus_h_eq)) {
          m_next_is_minus_previous = false;
          break;
        }
      }
    }
    return m_previous;
  }

  AsymIndex::AsymIndex(
    const sgtbx::SpaceGroup& SgOps,
    const sgtbx::ReciprocalSpaceASU& ASU,
    const Index& H)
  {
    m_TBF = SgOps.TBF();
    m_FriedelFlag = false;
    for(int iInv=0;iInv<SgOps.fInv();iInv++) {
      for(int iSMx=0;iSMx<SgOps.nSMx();iSMx++) {
        sgtbx::RTMx M = SgOps(0, iInv, iSMx);
        m_HR = H * M.Rpart();
        if (ASU.isInASU(m_HR)) {
          m_HT = sgtbx::HT_mod_1(H, M.Tpart());
          return;
        }
      }
    }
    cctbx_assert(!SgOps.isCentric());
    for(int iSMx=0;iSMx<SgOps.nSMx();iSMx++) {
      sgtbx::RTMx M = SgOps(0, 0, iSMx);
      m_HR = H * M.Rpart();
      if (ASU.isInASU(-m_HR)) {
        m_HT = sgtbx::HT_mod_1(H, M.Tpart());
        m_FriedelFlag = true;
        return;
      }
    }
    throw cctbx_internal_error();
  }

  AsymIndex::AsymIndex(SymEquivIndices const& SEMI)
  {
    m_TBF = SEMI[0].TBF();
    int iSelected = 0;
    Index SelectedH = SEMI[0].HR();
    m_FriedelFlag = false;
    for(int iList=0;iList<SEMI.N();iList++) {
      const SymEquivIndex& SEI = SEMI[iList];
      Index TrialH = SEI.HR();
      for(int iMate = 0; iMate < SEMI.fMates(true); iMate++) {
        if (iMate) TrialH = -TrialH;
        if (TrialH < SelectedH) {
          iSelected = iList;
          SelectedH = TrialH;
          m_FriedelFlag = (iMate != 0);
        }
      }
    }
    m_HR = SEMI[iSelected].HR();
    m_HT = SEMI[iSelected].HT();
  }

  AsymIndex::AsymIndex(
    const sgtbx::SpaceGroup& SgOps,
    const Index& H)
  {
    *this = AsymIndex(SymEquivIndices(SgOps, H));
  }

  join_sets::join_sets(
    af::shared<Index> miller_indices_0,
    af::shared<Index> miller_indices_1)
  : miller_indices_(miller_indices_0, miller_indices_1)
  {
    if (miller_indices_[0].id() == miller_indices_[1].id()) {
      // short-cut if same shared array
      pairs_.reserve(miller_indices_[0].size());
      for(std::size_t i=0;i<miller_indices_[0].size();i++) {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, i));
      }
      return;
    }
    typedef std::map<Index, std::size_t> lookup_map_type;
    lookup_map_type lookup_map;
    std::size_t i;
    for(i=0;i<miller_indices_[1].size();i++) {
      lookup_map[miller_indices_[1][i]] = i;
    }
    std::vector<bool> miller_indices_1_flags(miller_indices_[1].size(), false);
    for(i=0;i<miller_indices_[0].size();i++) {
      lookup_map_type::const_iterator
      l = lookup_map.find(miller_indices_[0][i]);
      if (l == lookup_map.end()) {
        singles_[0].push_back(i);
      }
      else {
        pairs_.push_back(af::tiny<std::size_t, 2>(i, l->second));
        miller_indices_1_flags[l->second] = true;
      }
    }
    for(i=0;i<miller_indices_[1].size();i++) {
      if (!miller_indices_1_flags[i]) singles_[1].push_back(i);
    }
  }

  af::shared<bool>
  join_sets::pair_selection(std::size_t i_array) const
  {
    size_assert_intrinsic();
    af::shared<bool> result(miller_indices_[i_array].size(), false);
    for(std::size_t i=0;i<pairs_.size();i++) {
      result[pairs_[i][i_array]] = true;
    }
    return result;
  }

  af::shared<bool>
  join_sets::single_selection(std::size_t i_array) const
  {
    size_assert_intrinsic();
    af::shared<bool> result(miller_indices_[i_array].size(), false);
    for(std::size_t i=0;i<singles_[i_array].size();i++) {
      result[singles_[i_array][i]] = true;
    }
    return result;
  }

  af::shared<Index>
  join_sets::paired_miller_indices(std::size_t i_array) const
  {
    size_assert_intrinsic();
    af::shared<Index> result;
    result.reserve(pairs_.size());
    for(std::size_t i=0;i<pairs_.size();i++) {
      result.push_back(miller_indices_[i_array][pairs_[i][i_array]]);
    }
    return result;
  }

  void join_bijvoet_mates::join_(sgtbx::ReciprocalSpaceASU const& asu)
  {
    typedef std::map<Index, std::size_t> lookup_map_type;
    lookup_map_type lookup_map;
    std::size_t i;
    for(i=0;i<miller_indices_.size();i++) {
      lookup_map[miller_indices_[i]] = i;
    }
    std::vector<bool> paired_already(miller_indices_.size(), false);
    for(i=0;i<miller_indices_.size();i++) {
      if (paired_already[i]) continue;
      lookup_map_type::const_iterator l = lookup_map.find(-miller_indices_[i]);
      if (l == lookup_map.end()) {
        singles_.push_back(i);
      }
      else {
        int asu_sign = asu.asu_sign(miller_indices_[i]);
        cctbx_assert(asu_sign != 0 || miller_indices_[i].is000());
        if (asu_sign > 0) {
          pairs_.push_back(af::tiny<std::size_t, 2>(i, l->second));
        }
        else {
          pairs_.push_back(af::tiny<std::size_t, 2>(l->second, i));
        }
        paired_already[l->second] = true;
      }
    }
  }

  af::shared<Index>
  join_bijvoet_mates::miller_indices_in_hemisphere(char plus_or_minus) const
  {
    cctbx_assert(plus_or_minus == '+' || plus_or_minus == '-');
    size_assert_intrinsic();
    std::size_t j = 0;
    if (plus_or_minus == '-') j = 1;
    af::shared<Index> result;
    result.reserve(pairs_.size());
    for(std::size_t i=0;i<pairs_.size();i++) {
      result.push_back(miller_indices_[pairs_[i][j]]);
    }
    return result;
  }

  binning::binning(
    uctbx::UnitCell const& unit_cell,
    std::size_t n_bins,
    af::shared<Index> miller_indices,
    double d_max,
    double d_min,
    double relative_tolerance)
  : unit_cell_(unit_cell)
  {
    if (!(d_max || d_min)) {
      af::double2 min_max_q = unit_cell.min_max_Q(miller_indices);
      if (!d_max && min_max_q[0]) d_max = 1 / std::sqrt(min_max_q[0]);
      if (!d_min && min_max_q[1]) d_min = 1 / std::sqrt(min_max_q[1]);
    }
    init_limits(n_bins, d_max, d_min, relative_tolerance);
  }

  void binning::init_limits(
    std::size_t n_bins,
    double d_max,
    double d_min,
    double relative_tolerance)
  {
    cctbx_assert(n_bins > 0);
    cctbx_assert(d_max >= 0);
    cctbx_assert(d_min > 0);
    cctbx_assert(d_min < d_max);
    double d_star_sq_min = 0;
    if (d_max) d_star_sq_min = 1 / (d_max * d_max);
    double     d_star_sq_max = 1 / (d_min * d_min);
    double span = d_star_sq_max - d_star_sq_min;
    d_star_sq_max += span * relative_tolerance;
    d_star_sq_min -= span * relative_tolerance;
    double r_low = std::sqrt(d_star_sq_min);
    double r_high = std::sqrt(d_star_sq_max);
    double volume_low = sphere_volume(r_low);
    double volume_per_bin = (sphere_volume(r_high) - volume_low) / n_bins;
    limits_.push_back(d_star_sq_min);
    for(std::size_t i_bin=1;i_bin<n_bins;i_bin++) {
      double r_sq_i = std::pow(
        (volume_low + i_bin * volume_per_bin) * 3 / constants::four_pi,
        2/3.);
      limits_.push_back(r_sq_i);
    }
    limits_.push_back(d_star_sq_max);
  }

  af::double2 binning::bin_d_range(std::size_t i_bin) const
  {
    return af::double2(bin_d_min(i_bin), bin_d_min(i_bin+1));
  }

  double binning::bin_d_min(std::size_t i_bin) const
  {
    if (i_bin == 0) return 0;
    if (i_bin == n_bins_all()) return 0;
    if (i_bin > n_bins_all()) throw error_index();
    return 1 / std::sqrt(limits_[i_bin - 1]);
  }

  std::size_t
  binning::get_i_bin(double d_star_sq) const
  {
    if (d_star_sq < limits_[0]) return 0;
    std::size_t i = 1;
    for(;i<limits_.size();i++) {
      if (d_star_sq < limits_[i]) return i;
    }
    return i;
  }

  binner::binner(binning const& bng, af::shared<Index> miller_indices)
  : binning(bng)
  {
    bin_indices_.reserve(miller_indices.size());
    for(std::size_t i=0;i<miller_indices.size();i++) {
      bin_indices_.push_back(this->get_i_bin(miller_indices[i]));
    }
  }

  std::size_t binner::count(std::size_t i_bin) const
  {
    cctbx_assert(i_bin < this->n_bins_all());
    std::size_t result = 0;
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      if (bin_indices_[i] == i_bin) result++;
    }
    return result;
  }

  af::shared<std::size_t> binner::counts() const
  {
    af::shared<std::size_t> result(this->n_bins_all());
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      std::size_t i_bin = bin_indices_[i];
      cctbx_assert(i_bin < result.size());
      result[i_bin]++;
    }
    return result;
  }

  af::shared<bool> binner::operator()(std::size_t i_bin) const
  {
    cctbx_assert(i_bin < this->n_bins_all());
    af::shared<bool> flags;
    flags.reserve(bin_indices_.size());
    for(std::size_t i=0;i<bin_indices_.size();i++) {
      flags.push_back(bin_indices_[i] == i_bin);
    }
    return flags;
  }

}} // namespace cctbx::miller
