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
#include <cctbx/sgtbx/groups.h>
#include <cctbx/sgtbx/miller_asu.h>

//#include <cctbx/array_family/simple_io.h> // XXX

namespace cctbx { namespace dmtbx {

  struct triplet_phase_relation
  {
    std::size_t weight;
    std::size_t ik;
    std::size_t ihmk;
    Miller::AsymIndex hmk;
  };

  template <typename FloatType>
  class triplet_invariants
  {
    private:
      typedef af::shared<triplet_phase_relation> tpr_array_type;

    public:
      triplet_invariants() {}

      triplet_invariants(sgtbx::SpaceGroupInfo const& SgInfo,
                         af::shared<Miller::Index> miller_indices,
                         af::shared<FloatType> e_values,
                         bool unique_triplet_phase_relations_only = false) // XXX not used
      {
        cctbx_assert(miller_indices.size() == e_values.size());
        Miller::index_span miller_index_span(miller_indices);
        typedef std::map<std::size_t, std::size_t> lookup_dict_type;
        lookup_dict_type lookup_dict;
        std::size_t i;
        for(i=0;i<miller_indices.size();i++) {
          lookup_dict[miller_index_span.pack(miller_indices[i])] = i;
        }
        cctbx_assert(lookup_dict.size() == miller_indices.size());
        sgtbx::ReciprocalSpaceASU asu(SgInfo);
        list_of_list_of_triplets_.reserve(miller_indices.size());
        triplet_phase_relation tpr;
        tpr.weight = 1;
        for(i=0;i<miller_indices.size();i++) {
          list_of_list_of_triplets_.push_back(tpr_array_type());
          Miller::Index h = miller_indices[i];
          for(tpr.ik=0;tpr.ik<miller_indices.size();tpr.ik++) {
            Miller::Index hmk = h - miller_indices[tpr.ik];
            tpr.hmk = Miller::AsymIndex(SgInfo.SgOps(), asu, hmk);
            Miller::Index asym_hmk = tpr.hmk.HermitianLayout().H();
            if (miller_index_span.is_in_domain(asym_hmk)) {
              typename lookup_dict_type::const_iterator
              ld_pos = lookup_dict.find(miller_index_span.pack(asym_hmk));
              if (ld_pos != lookup_dict.end()) {
                tpr.ihmk = ld_pos->second;
                cctbx_assert(miller_indices[tpr.ihmk] == asym_hmk);
                /*
                std::cout << h.ref()
                   << " " << miller_indices[tpr.ik].ref()
                   << " " << hmk.ref()
                   << " " << asym_hmk.ref()
                   << std::endl;
                */
                update_tpr_array(list_of_list_of_triplets_[i], tpr);
              }
            }
          }
        }
      }

      std::size_t total_number_of_triplets() const
      {
        std::size_t result = 0;
        for(std::size_t i=0;i<list_of_list_of_triplets_.size();i++) {
          result += list_of_list_of_triplets_[i].size();
        }
        return result;
      }

      FloatType average_number_of_triplets_per_reflection() const
      {
        return FloatType(total_number_of_triplets())
             / list_of_list_of_triplets_.size();
      }

      af::shared<FloatType>
      refine_phases(af::shared<Miller::Index> miller_indices,
                    af::shared<FloatType> e_values,
                    af::shared<FloatType> phases) const
      {
        cctbx_assert(miller_indices.size() == phases.size());
        cctbx_assert(miller_indices.size() == e_values.size());
        cctbx_assert(
          miller_indices.size() == list_of_list_of_triplets_.size());
        af::shared<FloatType> result;
        result.reserve(miller_indices.size());
        for(std::size_t i=0;i<list_of_list_of_triplets_.size();i++) {
          FloatType sum_sin(0);
          FloatType sum_cos(0);
          for(std::size_t j=0;j<list_of_list_of_triplets_[i].size();j++) {
            triplet_phase_relation const&
            tpr = list_of_list_of_triplets_[i][j];
            cctbx_assert(tpr.ik < miller_indices.size());
            cctbx_assert(tpr.ihmk < miller_indices.size());
            FloatType e_k = e_values[tpr.ik];
            FloatType phi_k = phases[tpr.ik];
            FloatType e_hmk = e_values[tpr.ihmk];
            FloatType phi_hmk = tpr.hmk.phase_in_rad(phases[tpr.ihmk]);
            FloatType w_e_k_e_hmk = tpr.weight * e_k * e_hmk;
            FloatType phi_k_phi_hmk = phi_k + phi_hmk;
            sum_sin += w_e_k_e_hmk * std::sin(phi_k_phi_hmk);
            sum_cos += w_e_k_e_hmk * std::cos(phi_k_phi_hmk);
          }
          cctbx_assert(sum_sin != 0 || sum_cos != 0);
          result.push_back(std::atan2(sum_sin, sum_cos));
        }
        return result;
      }

    private:
      af::shared<tpr_array_type> list_of_list_of_triplets_;

      static void update_tpr_array(
        tpr_array_type& tpr_array,
        triplet_phase_relation const& tpr)
      {
        for(std::size_t i=0;i<tpr_array.size();i++) {
          if (tpr_array[i].ik == tpr.ihmk && tpr_array[i].ihmk == tpr.ik) {
            tpr_array[i].weight++;
            return;
          }
        }
        tpr_array.push_back(tpr);
      }
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_TRIPLETS_H
