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
    std::size_t ik;
    std::size_t ihmk;
    Miller::AsymIndex asym_k;
    Miller::AsymIndex asym_hmk;

    triplet_phase_relation swap() const
    {
      triplet_phase_relation result;
      result.ik = ihmk;
      result.ihmk = ik;
      result.asym_k = asym_hmk;
      result.asym_hmk = asym_k;
      return result;
    }

    bool less(triplet_phase_relation const& other, bool eval_ht) const
    {
      return true;
      if (ik < other.ik) return true;
      if (ik > other.ik) return false;
      if (ihmk < other.ihmk) return true;
      if (ihmk > other.ihmk) return false;
      if (!eval_ht) return false;
      if (asym_k.HT() < other.asym_k.HT()) return true;
      if (asym_k.HT() > other.asym_k.HT()) return false;
      if (!asym_k.FriedelFlag() && other.asym_k.FriedelFlag()) return true;
      if (asym_k.FriedelFlag() && !other.asym_k.FriedelFlag()) return false;
      if (asym_hmk.HT() < other.asym_hmk.HT()) return true;
      if (asym_hmk.HT() > other.asym_hmk.HT()) return false;
      if (!asym_hmk.FriedelFlag() && other.asym_hmk.FriedelFlag()) return true;
      return false;
    }

    bool operator<(triplet_phase_relation const& other) const
    {
      return less(other, true);
    }

    struct less_no_eval_ht
    {
      bool operator()(
        triplet_phase_relation const& lhs,
        triplet_phase_relation const& rhs)
      {
        return lhs.less(rhs, false);
      }
    };
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
      {
        cctbx_assert(miller_indices.size() == e_values.size());
        sgtbx::ReciprocalSpaceASU asu(SgInfo);
        std::size_t i=0;
        // Assert that all Miller indices are in the standard asymmetric unit.
        for(i=0;i<miller_indices.size();i++) {
          cctbx_assert(
            Miller::AsymIndex(SgInfo.SgOps(), asu, miller_indices[i]
              ).one_column(true).H()
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
        triplet_phase_relation tpr;
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          //if (ih == tpr.ik) continue;
          Miller::Index h = miller_indices[ih];
          for(tpr.ik=0;tpr.ik<miller_indices.size();tpr.ik++) {
            Miller::Index k = miller_indices[tpr.ik];
            sgtbx::SymEquivMillerIndices
            sym_eq_k = SgInfo.SgOps().getEquivMillerIndices(k);
            for (std::size_t ik_eq=0;ik_eq<sym_eq_k.M(true);ik_eq++) {
              Miller::Index k_eq = sym_eq_k(ik_eq).H();
              Miller::Index hmk = h - k_eq;
              //std::cout << h.ref()
              //   << " " << k_eq.ref()
              //   << " " << hmk.ref()
              //   << " under consideration"
              //   << std::endl;
              tpr.asym_hmk = Miller::AsymIndex(SgInfo.SgOps(), asu, hmk);
              Miller::Index asym_hmk = tpr.asym_hmk.one_column(true).H();
              //std::cout << asym_hmk.ref()
              //   << " asym_hmk"
              //   << std::endl;
              if (miller_index_span.is_in_domain(asym_hmk)) {
                typename lookup_dict_type::const_iterator
                ld_pos = lookup_dict.find(miller_index_span.pack(asym_hmk));
                if (ld_pos != lookup_dict.end()) {
                  tpr.ihmk = ld_pos->second;
                  cctbx_assert(miller_indices[tpr.ihmk] == asym_hmk);
                  //if (tpr.ihmk == ih) continue;
                  //if (tpr.ihmk == tpr.ik) continue;
                  tpr.asym_k = Miller::AsymIndex(SgInfo.SgOps(), asu, k_eq);
                  //std::cout << h.ref()
                  //   << " " << k_eq.ref()
                  //   << " " << hmk.ref()
                  //   << std::endl;
                  if (tpr.ik > tpr.ihmk) {
                    //list_of_tpr_maps_[ih][tpr.swap()]++;
                    list_of_tpr_maps_[ih][tpr]++;
                  }
                  else {
                    list_of_tpr_maps_[ih][tpr]++;
                  }
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
               << " " << miller_indices[tpr.ik].ref()
               << " " << miller_indices[tpr.ihmk].ref()
               << " " << lij->second
               << std::endl;
          }
        }
      }

      triplet_invariants
      unique_triplets()
      {
        typedef
          std::map<
            triplet_phase_relation,
            std::size_t,
            triplet_phase_relation::less_no_eval_ht> umap_type;
        triplet_invariants result;
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        for(std::size_t i=0;i<list_of_tpr_maps_.size();i++,li++) {
          umap_type umap;
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            umap[lij->first] += lij->second;
          }
          tpr_map_type new_map;
          for(umap_type::const_iterator ui=umap.begin();ui!=umap.end();ui++) {
            new_map[ui->first] = ui->second;
          }
          result.list_of_tpr_maps_.append(new_map);
        }
        return result;
      }

      af::shared<FloatType>
      apply_tangent_formula(af::shared<FloatType> e_values,
                            af::shared<FloatType> phases,
                            bool ignore_weights = false) const
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
            cctbx_assert(tpr.ik < e_values.size());
            cctbx_assert(tpr.ihmk < e_values.size());
            FloatType e_k = e_values[tpr.ik];
            FloatType phi_k = tpr.asym_k.phase_in(phases[tpr.ik]);
            FloatType e_hmk = e_values[tpr.ihmk];
            FloatType phi_hmk = tpr.asym_hmk.phase_in(phases[tpr.ihmk]);
            FloatType e_k_e_hmk = e_k * e_hmk;
            if (!ignore_weights) e_k_e_hmk *= lij->second;
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
                      af::shared<FloatType> phases,
                      bool ignore_weights = false) const
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
            cctbx_assert(tpr.ik < e_values.size());
            cctbx_assert(tpr.ihmk < e_values.size());
            FloatType e_k = e_values[tpr.ik];
            FloatType phi_k = tpr.asym_k.phase_in(phases[tpr.ik]);
            FloatType e_hmk = e_values[tpr.ihmk];
            FloatType phi_hmk = tpr.asym_hmk.phase_in(phases[tpr.ihmk]);
            std::complex<FloatType> e_k_complex = std::polar(e_k, phi_k);
            std::complex<FloatType> e_hmk_complex = std::polar(e_hmk, phi_hmk);
            if (!ignore_weights) {
              estimated_e_h
                += FloatType(lij->second) * (e_k_complex * e_hmk_complex);
            }
            else {
              estimated_e_h
                +=                          (e_k_complex * e_hmk_complex);
            }
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
      list_of_tpr_maps_type list_of_tpr_maps_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_TRIPLETS_H
