/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created based on triplet.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_DMTBX_TRIPLET_PHASE_RELATION_H
#define CCTBX_DMTBX_TRIPLET_PHASE_RELATION_H

#include <cctbx/math/mod.h>
#include <cctbx/error.h>
#include <scitbx/constants.h>

namespace cctbx { namespace dmtbx {

  class triplet_phase_relation
  {
    public:
      triplet_phase_relation() {}

      //! Not available in Python.
      triplet_phase_relation(
        std::size_t ik,
        bool friedel_flag_k,
        int ht_k,
        std::size_t ihmk,
        bool friedel_flag_hmk,
        int ht_hmk,
        int t_den)
      {
        CCTBX_ASSERT(ht_k >= 0);
        CCTBX_ASSERT(ht_hmk >= 0);
        if (ik < ihmk || (ik == ihmk && friedel_flag_k == false)) {
          ik_ = ik;
          ihmk_ = ihmk;
          friedel_flag_k_ = friedel_flag_k;
          friedel_flag_hmk_ = friedel_flag_hmk;
        }
        else {
          ik_ = ihmk;
          ihmk_ = ik;
          friedel_flag_k_ = friedel_flag_hmk;
          friedel_flag_hmk_ = friedel_flag_k;
        }
        if (!friedel_flag_k) ht_k *= -1;
        if (!friedel_flag_hmk) ht_hmk *= -1;
        ht_sum_ = math::mod_positive(ht_k + ht_hmk, t_den);
      }

      std::size_t
      ik() const { return ik_; }

      bool
      friedel_flag_k() const { return friedel_flag_k_; }

      std::size_t
      ihmk() const { return ihmk_; }

      bool
      friedel_flag_hmk() const { return friedel_flag_hmk_; }

      int
      ht_sum() const { return ht_sum_; }

      //! Not available in Python.
      bool operator<(triplet_phase_relation const& other) const
      {
        if (ik_ < other.ik_) return true;
        if (ik_ > other.ik_) return false;
        if (ihmk_ < other.ihmk_) return true;
        if (ihmk_ > other.ihmk_) return false;
        if (ht_sum_ < other.ht_sum_) return true;
        if (ht_sum_ > other.ht_sum_) return false;
        if (!friedel_flag_k_ && other.friedel_flag_k_) return true;
        if (friedel_flag_k_ && !other.friedel_flag_k_) return false;
        if (!friedel_flag_hmk_ && other.friedel_flag_hmk_) return true;
        return false;
      }

      //! Not available in Python.
      template <typename FloatType>
      FloatType
      phi_k_phi_hmk(const FloatType* phases, int t_den) const
      {
        FloatType phi_k = phases[ik_];
        if (friedel_flag_k_) phi_k = -phi_k;
        FloatType phi_hmk = phases[ihmk_];
        if (friedel_flag_hmk_) phi_hmk = -phi_hmk;
        return phi_k + phi_hmk + (scitbx::constants::two_pi * ht_sum_) / t_den;
      }

    protected:
      std::size_t ik_;
      bool friedel_flag_k_;
      std::size_t ihmk_;
      bool friedel_flag_hmk_;
      int ht_sum_;
  };

  class weighted_triplet_phase_relation : public triplet_phase_relation
  {
    public:
      weighted_triplet_phase_relation() {}

      weighted_triplet_phase_relation(
        triplet_phase_relation const& tpr,
        std::size_t weight)
      :
        triplet_phase_relation(tpr),
        weight_(weight)
      {}

      std::size_t
      weight() const { return weight_; }

    protected:
      std::size_t weight_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_TRIPLET_PHASE_RELATION_H
