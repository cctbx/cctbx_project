#ifndef CCTBX_DMTBX_TRIPLET_PHASE_RELATION_H
#define CCTBX_DMTBX_TRIPLET_PHASE_RELATION_H

#include <scitbx/math/modulo.h>
#include <cctbx/error.h>
#include <scitbx/constants.h>

namespace cctbx { namespace dmtbx {

  //! Triplet phase relation h = k + h-k.
  class triplet_phase_relation
  {
    public:
      //! Default constructor. Some data members are not initialized!
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
        ht_sum_ = scitbx::math::mod_positive(ht_k + ht_hmk, t_den);
      }

      //! Pointer to unique Miller index for k.
      std::size_t
      ik() const { return ik_; }

      //! Friedel flag for k.
      bool
      friedel_flag_k() const { return friedel_flag_k_; }

      //! Pointer to unique Miller index for h-k.
      std::size_t
      ihmk() const { return ihmk_; }

      //! Friedel flag for h-k.
      bool
      friedel_flag_hmk() const { return friedel_flag_hmk_; }

      //! Sum of the phase shifts due non-zero translation parts.
      int
      ht_sum() const { return ht_sum_; }

      //! True if sigma-2 relation.
      /*! A triplet phase relation h = k + h-k is a "sigma-2" relation
          only if the three Miller indices involved are not related
          by symmetry.

          ih is a pointer to the unique Miller index for h.
       */
      bool
      is_sigma_2(std::size_t ih) const
      {
        return ik_ != ih && ihmk_ != ih && ik_ != ihmk_;
      }

      //! True if ik()==other.ik() and ihmk()==other.ihmk().
      bool is_similar_to(triplet_phase_relation const& other) const
      {
        return ik_ == other.ik_ && ihmk_ == other.ihmk_;
      }

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

  //! Weighted triplet phase relation.
  class weighted_triplet_phase_relation : public triplet_phase_relation
  {
    public:
      //! Default constructor. Some data members are not initialized!
      weighted_triplet_phase_relation() {}

      //! Not available in Python.
      weighted_triplet_phase_relation(
        triplet_phase_relation const& tpr,
        std::size_t weight)
      :
        triplet_phase_relation(tpr),
        weight_(weight)
      {}

      //! Number of times the triplet occurs.
      std::size_t
      weight() const { return weight_; }

    protected:
      std::size_t weight_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_TRIPLET_PHASE_RELATION_H
