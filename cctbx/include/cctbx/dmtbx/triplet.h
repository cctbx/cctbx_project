/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

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

#ifndef CCTBX_DMTBX_TRIPLET_H
#define CCTBX_DMTBX_TRIPLET_H

#include <cctbx/miller/asu.h>
#include <cctbx/miller/index_span.h>
#include <cctbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <map>

namespace cctbx { namespace dmtbx {

  struct triplet_phase_relation
  {
    triplet_phase_relation(
      std::size_t ik,
      std::size_t ihmk,
      miller::sym_equiv_index const& k_seq,
      miller::asym_index const& asym_hmk)
    {
      CCTBX_ASSERT(k_seq.t_den() == asym_hmk.t_den());
      CCTBX_ASSERT(k_seq.ht() >= 0);
      CCTBX_ASSERT(asym_hmk.ht() >= 0);
      if (ik <= ihmk) {
        ik_ = ik;
        ihmk_ = ihmk;
        friedel_flag_k_ = k_seq.friedel_flag();
        friedel_flag_hmk_ = asym_hmk.friedel_flag();
      }
      else {
        ik_ = ihmk;
        ihmk_ = ik;
        friedel_flag_k_ = asym_hmk.friedel_flag();
        friedel_flag_hmk_ = k_seq.friedel_flag();
      }
      int k_ht = k_seq.ht();
      if (k_seq.friedel_flag()) k_ht *= -1;
      ht_sum_ = math::mod_positive(-k_ht + asym_hmk.ht(), k_seq.t_den());
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
    phi_k_phi_hmk(const FloatType* phases, int t_den) const
    {
      FloatType phi_k = phases[ik_];
      if (friedel_flag_k_) phi_k = -phi_k;
      FloatType phi_hmk = phases[ihmk_];
      if (friedel_flag_hmk_) phi_hmk = -phi_hmk;
      return phi_k + phi_hmk + (scitbx::constants::two_pi * ht_sum_) / t_den;
    }

    std::size_t ik_;
    std::size_t ihmk_;
    bool friedel_flag_k_;
    bool friedel_flag_hmk_;
    int ht_sum_;
  };

  template <typename FloatType = double>
  class triplet_invariants
  {
    private:
      typedef std::map<triplet_phase_relation, std::size_t> tpr_map_type;
      typedef af::shared<tpr_map_type> list_of_tpr_maps_type;

    public:
      triplet_invariants() {}

      triplet_invariants(
        sgtbx::space_group_type const& sg_type,
        af::const_ref<miller::index<> > const& miller_indices,
        bool sigma_2,
        bool other_than_sigma_2)
      :
        t_den_(sg_type.group().t_den())
      {
        sgtbx::reciprocal_space::asu asu(sg_type);
        std::size_t i=0;
        // Assert that all Miller indices are in the standard asymmetric unit.
        for(i=0;i<miller_indices.size();i++) {
          CCTBX_ASSERT(
               miller::asym_index(sg_type.group(), asu, miller_indices[i]).h()
            == miller_indices[i]);
        }
        miller::index_span miller_index_span(miller_indices);
        typedef std::map<std::size_t, std::size_t> lookup_dict_type;
        lookup_dict_type lookup_dict;
        for(i=0;i<miller_indices.size();i++) {
          lookup_dict[miller_index_span.pack(miller_indices[i])] = i;
        }
        CCTBX_ASSERT(lookup_dict.size() == miller_indices.size());
        list_of_tpr_maps_.reserve(miller_indices.size());
        for(i=0;i<miller_indices.size();i++) {
          list_of_tpr_maps_.push_back(tpr_map_type());
        }
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          miller::index<> h = miller_indices[ih];
          for(std::size_t ik=0;ik<miller_indices.size();ik++) {
            if (ik == ih) { if (!other_than_sigma_2) continue; }
            else          { if (!sigma_2) continue; }
            miller::index<> k = miller_indices[ik];
            miller::sym_equiv_indices sym_eq_k(sg_type.group(), k);
            int mult = sym_eq_k.multiplicity(false);
            for(std::size_t ik_eq=0;ik_eq<mult;ik_eq++) {
              miller::sym_equiv_index k_seq = sym_eq_k(ik_eq);
              miller::index<> k_eq = k_seq.h();
              miller::index<> hmk = h - k_eq;
              miller::asym_index asym_hmk(sg_type.group(), asu, hmk);
              miller::index<> asym_hmk_h = asym_hmk.h();
              if (miller_index_span.is_in_domain(asym_hmk_h)) {
                typename lookup_dict_type::const_iterator
                ld_pos = lookup_dict.find(miller_index_span.pack(asym_hmk_h));
                if (ld_pos != lookup_dict.end()) {
                  std::size_t ihmk = ld_pos->second;
                  CCTBX_ASSERT(miller_indices[ihmk] == asym_hmk_h);
                  if (ihmk == ih || ihmk == ik) {
                    if (!other_than_sigma_2) continue;
                  }
                  else {
                    if (!sigma_2) continue;
                  }
                  triplet_phase_relation tpr(ik, ihmk, k_seq, asym_hmk);
                  list_of_tpr_maps_[ih][tpr]++;
                }
              }
            }
          }
        }
      }

      std::size_t
      number_of_weighted_triplets() const
      {
        std::size_t result = 0;
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++) {
          result += list_of_tpr_maps_[i].size();
        }
        return result;
      }

      std::size_t
      total_number_of_triplets() const
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

      FloatType
      average_number_of_triplets_per_reflection() const
      {
        return FloatType(total_number_of_triplets())
             / list_of_tpr_maps_.size();
      }

      void
      dump_triplets(af::const_ref<miller::index<> > const& miller_indices)
      {
        CCTBX_ASSERT(miller_indices.size() == list_of_tpr_maps_.size());
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++,li++) {
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            std::cout << miller_indices[i].const_ref()
               << " " << miller_indices[tpr.ik_].const_ref()
               << " " << tpr.friedel_flag_k_
               << " " << miller_indices[tpr.ihmk_].const_ref()
               << " " << tpr.friedel_flag_hmk_
               << " shift: " << tpr.ht_sum_
               << " w: " << lij->second
               << std::endl;
          }
        }
      }

      af::shared<std::size_t>
      n_relations(bool discard_weights,
                  bool first_only) const
      {
        af::shared<std::size_t>
          result((af::reserve(list_of_tpr_maps_.size())));
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        std::size_t n_miller_indices = list_of_tpr_maps_.size();
        std::size_t prev_ik = n_miller_indices;
        std::size_t prev_ihmk = n_miller_indices;
        for(std::size_t i_h=0;i_h<n_miller_indices;i_h++,li++) {
          std::size_t n = 0;
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            if (first_only) {
              if (tpr.ik_ == prev_ik && tpr.ihmk_ == prev_ihmk) continue;
              prev_ik = tpr.ik_;
              prev_ihmk = tpr.ihmk_;
            }
            if (!discard_weights) n += lij->second;
            else n++;
          }
          result.push_back(n);
        }
        return result;
      }

      af::shared<FloatType>
      sum_of_e_products(
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<FloatType> const& e_values,
        bool discard_weights,
        bool first_only) const
      {
        CCTBX_ASSERT(miller_indices.size() == list_of_tpr_maps_.size());
        CCTBX_ASSERT(miller_indices.size() == e_values.size());
        af::shared<FloatType> result;
        result.reserve(e_values.size());
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        std::size_t n_miller_indices = list_of_tpr_maps_.size();
        std::size_t prev_ik = n_miller_indices;
        std::size_t prev_ihmk = n_miller_indices;
        for(std::size_t i_h=0;i_h<n_miller_indices;i_h++,li++) {
          FloatType sum = 0;
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            if (first_only) {
              if (tpr.ik_ == prev_ik && tpr.ihmk_ == prev_ihmk) continue;
              prev_ik = tpr.ik_;
              prev_ihmk = tpr.ihmk_;
            }
            FloatType e_k = e_values[tpr.ik_];
            FloatType e_hmk = e_values[tpr.ihmk_];
            FloatType ee = e_k * e_hmk;
            if (!discard_weights) ee *= lij->second;
            sum += ee;
          }
          result.push_back(sum);
        }
        return result;
      }

      af::shared<FloatType>
      apply_tangent_formula(
        af::const_ref<FloatType> const& e_values,
        af::const_ref<FloatType> const& phases,
        af::const_ref<bool> const& selection_fixed,
        af::const_ref<std::size_t> const& extrapolation_order,
        bool reuse_results,
        bool discard_weights,
        bool first_only) const
      {
        FloatType sum_cutoff(1.e-10); // XXX
        CCTBX_ASSERT(e_values.size() == list_of_tpr_maps_.size());
        CCTBX_ASSERT(phases.size() == e_values.size());
        CCTBX_ASSERT(selection_fixed.size() == e_values.size());
        CCTBX_ASSERT(   extrapolation_order.size() == 0
                     || extrapolation_order.size() == e_values.size());
        CCTBX_ASSERT(first_only == false || discard_weights == true);
        af::shared<FloatType> result;
        result.assign(phases.begin(), phases.end());
        const FloatType* phase_source = (
          reuse_results ? result.begin() : phases.begin());
        std::vector<bool> fixed_or_extrapolated(
          selection_fixed.begin(), selection_fixed.end());
        list_of_tpr_maps_type::const_iterator
          li_begin = list_of_tpr_maps_.begin();
        std::size_t prev_ik = e_values.size();
        std::size_t prev_ihmk = e_values.size();
        std::size_t i_h;
        for(std::size_t i_p=0;i_p<phases.size();i_p++) {
          if (extrapolation_order.size() == 0) {
            i_h = i_p;
          }
          else {
            i_h = extrapolation_order[i_p];
            CCTBX_ASSERT(i_h < e_values.size());
          }
          if (selection_fixed[i_h]) continue;
          CCTBX_ASSERT(!fixed_or_extrapolated[i_h]);
          list_of_tpr_maps_type::const_iterator li = li_begin + i_h;
          FloatType sum_sin(0);
          FloatType sum_cos(0);
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            if (first_only) {
              if (tpr.ik_ == prev_ik && tpr.ihmk_ == prev_ihmk) continue;
              prev_ik = tpr.ik_;
              prev_ihmk = tpr.ihmk_;
            }
            CCTBX_ASSERT(tpr.ik_ < e_values.size());
            CCTBX_ASSERT(tpr.ihmk_ < e_values.size());
            if (reuse_results) {
              if (!fixed_or_extrapolated[tpr.ik_]) continue;
              if (!fixed_or_extrapolated[tpr.ihmk_]) continue;
            }
            FloatType e_k = e_values[tpr.ik_];
            FloatType e_hmk = e_values[tpr.ihmk_];
            FloatType e_k_e_hmk = e_k * e_hmk;
            if (!discard_weights) e_k_e_hmk *= lij->second;
            FloatType phi_k_phi_hmk = tpr.phi_k_phi_hmk(phase_source, t_den_);
            sum_sin += e_k_e_hmk * std::sin(phi_k_phi_hmk);
            sum_cos += e_k_e_hmk * std::cos(phi_k_phi_hmk);
          }
          if (   scitbx::fn::absolute(sum_sin) >= sum_cutoff
              || scitbx::fn::absolute(sum_cos) >= sum_cutoff) {
            result[i_h] = std::atan2(sum_sin, sum_cos);
            fixed_or_extrapolated[i_h] = true;
          }
        }
        return result;
      }

    private:
      int t_den_;
      list_of_tpr_maps_type list_of_tpr_maps_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_TRIPLET_H
