// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created (R.W. Grosse-Kunstleve)
 */

/* Reference:
     C.M. Weeks, P.D. Adams, J. Berendzen, A.T. Brunger, E.J. Dodson,
     R.W. Grosse-Kunstleve, T.R. Schneider, G.M. Sheldrick,
     T.C. Terwilliger, M. Turkenburg, I. Uson
     Automatic solution of heavy-atom substructures.
     Methods in Enzymology, in press.
 */

#ifndef CCTBX_TRIPLETS_H
#define CCTBX_TRIPLETS_H

#include <map>
#include <cctbx/error.h>
#include <cctbx/array_family/shared.h>
#include <cctbx/array_family/simple_io.h>
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_asu.h>

namespace cctbx { namespace dmtbx {

  struct triplet_phase_relation
  {
    triplet_phase_relation(
      std::size_t ik,
      std::size_t ihmk,
      Miller::SymEquivIndex const& k_seq,
      Miller::AsymIndex const& asym_hmk)
    {
      cctbx_assert(k_seq.TBF() == asym_hmk.TBF());
      cctbx_assert(k_seq.HT() >= 0);
      cctbx_assert(asym_hmk.HT() >= 0);
      if (ik <= ihmk) {
        ik_ = ik;
        ihmk_ = ihmk;
        friedel_flag_k_ = k_seq.FriedelFlag();
        friedel_flag_hmk_ = asym_hmk.FriedelFlag();
      }
      else {
        ik_ = ihmk;
        ihmk_ = ik;
        friedel_flag_k_ = asym_hmk.FriedelFlag();
        friedel_flag_hmk_ = k_seq.FriedelFlag();
      }
      int k_ht = k_seq.HT();
      if (k_seq.FriedelFlag()) k_ht *= -1;
      ht_sum_ = sgtbx::modPositive(-k_ht + asym_hmk.HT(), k_seq.TBF());
    }

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

    template <typename FloatType>
    FloatType
    phi_k_phi_hmk(const FloatType* phases, int TBF) const
    {
      FloatType phi_k = phases[ik_];
      if (friedel_flag_k_) phi_k = -phi_k;
      FloatType phi_hmk = phases[ihmk_];
      if (friedel_flag_hmk_) phi_hmk = -phi_hmk;
      return phi_k + phi_hmk + (constants::two_pi * ht_sum_) / TBF;
    }

    std::size_t ik_;
    std::size_t ihmk_;
    bool friedel_flag_k_;
    bool friedel_flag_hmk_;
    int ht_sum_;
  };

  template <typename FloatType>
  class triplet_invariants
  {
    private:
      typedef std::map<triplet_phase_relation, std::size_t> tpr_map_type;
      typedef af::shared<tpr_map_type> list_of_tpr_maps_type;

    public:
      triplet_invariants() {}

      triplet_invariants(sgtbx::SpaceGroupInfo const& SgInfo,
                         af::shared<Miller::Index> miller_indices,
                         af::shared<FloatType> e_values)
      : TBF_(SgInfo.SgOps().TBF())
      {
        cctbx_assert(miller_indices.size() == e_values.size());
        sgtbx::ReciprocalSpaceASU asu(SgInfo);
        std::size_t i=0;
        // Assert that all Miller indices are in the standard asymmetric unit.
        for(i=0;i<miller_indices.size();i++) {
          cctbx_assert(
               Miller::AsymIndex(SgInfo.SgOps(), asu, miller_indices[i]).H()
            == miller_indices[i]);
        }
        Miller::index_span miller_index_span(miller_indices);
        typedef std::map<std::size_t, std::size_t> lookup_dict_type;
        lookup_dict_type lookup_dict;
        for(i=0;i<miller_indices.size();i++) {
          lookup_dict[miller_index_span.pack(miller_indices[i])] = i;
        }
        cctbx_assert(lookup_dict.size() == miller_indices.size());
        list_of_tpr_maps_.reserve(miller_indices.size());
        for(i=0;i<miller_indices.size();i++) {
          list_of_tpr_maps_.push_back(tpr_map_type());
        }
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          //if (ik == ih) continue;
          Miller::Index h = miller_indices[ih];
          for(std::size_t ik=0;ik<miller_indices.size();ik++) {
            Miller::Index k = miller_indices[ik];
            sgtbx::SymEquivMillerIndices
            sym_eq_k = SgInfo.SgOps().getEquivMillerIndices(k);
            for (std::size_t ik_eq=0;ik_eq<sym_eq_k.M(true);ik_eq++) {
              Miller::SymEquivIndex k_seq = sym_eq_k(ik_eq);
              Miller::Index k_eq = k_seq.H();
              Miller::Index hmk = h - k_eq;
              Miller::AsymIndex asym_hmk(SgInfo.SgOps(), asu, hmk);
              Miller::Index asym_hmk_h = asym_hmk.H();
              if (miller_index_span.is_in_domain(asym_hmk_h)) {
                typename lookup_dict_type::const_iterator
                ld_pos = lookup_dict.find(miller_index_span.pack(asym_hmk_h));
                if (ld_pos != lookup_dict.end()) {
                  std::size_t ihmk = ld_pos->second;
                  cctbx_assert(miller_indices[ihmk] == asym_hmk_h);
                  //if (ihmk == ih) continue;
                  //if (ihmk == ik) continue;
                  triplet_phase_relation tpr(ik, ihmk, k_seq, asym_hmk);
                  list_of_tpr_maps_[ih][tpr]++;
                }
              }
            }
          }
        }
      }

      std::size_t number_of_weighted_triplets() const
      {
        std::size_t result = 0;
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++) {
          result += list_of_tpr_maps_[i].size();
        }
        return result;
      }

      std::size_t total_number_of_triplets() const
      {
        std::size_t result = 0;
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++,li++) {
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            result += lij->second;
          }
        }
        return result;
      }

      FloatType average_number_of_triplets_per_reflection() const
      {
        return FloatType(total_number_of_triplets())
             / list_of_tpr_maps_.size();
      }

      std::size_t n_relations(std::size_t ih) const
      {
        std::size_t result = 0;
        cctbx_assert(ih < list_of_tpr_maps_.size());
        list_of_tpr_maps_type::const_iterator
        li = list_of_tpr_maps_.begin() + ih;
        for(tpr_map_type::const_iterator
            lij=li->begin();lij!=li->end();lij++) {
          result += lij->second;
        }
        return result;
      }

      void dump_triplets(af::shared<Miller::Index> miller_indices)
      {
        cctbx_assert(
          miller_indices.size() == list_of_tpr_maps_.size());
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++,li++) {
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            std::cout << miller_indices[i].ref()
               << " " << miller_indices[tpr.ik_].ref()
               << " " << miller_indices[tpr.ihmk_].ref()
               << " " << lij->second
               << std::endl;
          }
        }
      }

      af::shared<FloatType>
      apply_tangent_formula(af::shared<FloatType> e_values,
                            af::shared<FloatType> phases) const
      {
        FloatType sum_cutoff(1.e-10); // XXX
        cctbx_assert(e_values.size() == list_of_tpr_maps_.size());
        cctbx_assert(e_values.size() == phases.size());
        af::shared<FloatType> result;
        result.reserve(e_values.size());
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++,li++) {
          FloatType sum_sin(0);
          FloatType sum_cos(0);
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const&
            tpr = lij->first;
            cctbx_assert(tpr.ik_ < e_values.size());
            cctbx_assert(tpr.ihmk_ < e_values.size());
            FloatType e_k = e_values[tpr.ik_];
            FloatType e_hmk = e_values[tpr.ihmk_];
            FloatType e_k_e_hmk = lij->second * e_k * e_hmk;
            FloatType phi_k_phi_hmk = tpr.phi_k_phi_hmk(phases.begin(), TBF_);
            sum_sin += e_k_e_hmk * std::sin(phi_k_phi_hmk);
            sum_cos += e_k_e_hmk * std::cos(phi_k_phi_hmk);
          }
          if (   math::abs(sum_sin) < sum_cutoff
              && math::abs(sum_cos) < sum_cutoff) {
            result.push_back(phases[i]);
          }
          else {
            result.push_back(std::atan2(sum_sin, sum_cos));
          }
        }
        return result;
      }

    private:
      int TBF_;
      list_of_tpr_maps_type list_of_tpr_maps_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_TRIPLETS_H
