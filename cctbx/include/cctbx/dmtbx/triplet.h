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
      Miller::AsymIndex asym_k,
      Miller::AsymIndex asym_hmk)
    {
      if (ik <= ihmk) {
        ik_ = ik;
        ihmk_ = ihmk;
        asym_k_ = asym_k;
        asym_hmk_ = asym_hmk;
      }
      else {
        ik_ = ihmk;
        ihmk_ = ik;
        asym_k_ = asym_hmk;
        asym_hmk_ = asym_k;
      }
    }

    bool operator<(triplet_phase_relation const& other) const
    {
      if (ik_ < other.ik_) return true;
      if (ik_ > other.ik_) return false;
      if (ihmk_ < other.ihmk_) return true;
      if (ihmk_ > other.ihmk_) return false;
      if (asym_k_.HT() < other.asym_k_.HT()) return true;
      if (asym_k_.HT() > other.asym_k_.HT()) return false;
      if (!asym_k_.FriedelFlag() && other.asym_k_.FriedelFlag()) return true;
      if (asym_k_.FriedelFlag() && !other.asym_k_.FriedelFlag()) return false;
      if (asym_hmk_.HT() < other.asym_hmk_.HT()) return true;
      if (asym_hmk_.HT() > other.asym_hmk_.HT()) return false;
      if (!asym_hmk_.FriedelFlag() && other.asym_hmk_.FriedelFlag())return true;
      return false;
    }

    std::size_t ik_;
    std::size_t ihmk_;
    Miller::AsymIndex asym_k_;
    Miller::AsymIndex asym_hmk_;
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
              Miller::Index k_eq = sym_eq_k(ik_eq).H();
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
                  Miller::AsymIndex asym_k(SgInfo.SgOps(), asu, k_eq);
                  triplet_phase_relation tpr(ik, ihmk, asym_k, asym_hmk);
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
            FloatType phi_k = tpr.asym_k_.phase_in(phases[tpr.ik_]);
            FloatType e_hmk = e_values[tpr.ihmk_];
            FloatType phi_hmk = tpr.asym_hmk_.phase_in(phases[tpr.ihmk_]);
            FloatType e_k_e_hmk = lij->second * e_k * e_hmk;
            FloatType phi_k_phi_hmk = phi_k + phi_hmk;
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

      af::shared<FloatType>
      estimate_phases(af::shared<FloatType> e_values,
                      af::shared<FloatType> phases) const
      {
        FloatType estimated_e_h_cutoff(1.e-10); // XXX
        cctbx_assert(e_values.size() == list_of_tpr_maps_.size());
        cctbx_assert(e_values.size() == phases.size());
        af::shared<FloatType> result;
        result.reserve(e_values.size());
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++,li++) {
          std::complex<FloatType> estimated_e_h(0);
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const&
            tpr = lij->first;
            cctbx_assert(tpr.ik_ < e_values.size());
            cctbx_assert(tpr.ihmk_ < e_values.size());
            FloatType e_k = e_values[tpr.ik_];
            FloatType phi_k = tpr.asym_k_.phase_in(phases[tpr.ik_]);
            FloatType e_hmk = e_values[tpr.ihmk_];
            FloatType phi_hmk = tpr.asym_hmk_.phase_in(phases[tpr.ihmk_]);
            std::complex<FloatType> e_k_complex = std::polar(e_k, phi_k);
            std::complex<FloatType> e_hmk_complex = std::polar(e_hmk, phi_hmk);
            estimated_e_h
              += FloatType(lij->second) * (e_k_complex * e_hmk_complex);
          }
          if (std::abs(estimated_e_h) < estimated_e_h_cutoff) {
            result.push_back(phases[i]);
          }
          else {
            result.push_back(std::arg(estimated_e_h));
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
