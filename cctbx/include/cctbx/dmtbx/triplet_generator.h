/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created based on triplet.h (R.W. Grosse-Kunstleve)
 */

/* References:
     C.M. Weeks, P.D. Adams, J. Berendzen, A.T. Brunger, E.J. Dodson,
     R.W. Grosse-Kunstleve, T.R. Schneider, G.M. Sheldrick,
     T.C. Terwilliger, M. Turkenburg, I. Uson
     Automatic solution of heavy-atom substructures.
     Methods in Enzymology, in press.

     Sheldrick, G.M. (1982).
     Crystallographic algorithms for mini- and maxi-computers.
     In: Computational Crystallography, Ed. D. Sayre,
     Oxford University Press, 506-514.
 */

#ifndef CCTBX_DMTBX_TRIPLET_GENERATOR_H
#define CCTBX_DMTBX_TRIPLET_GENERATOR_H

#include <cctbx/dmtbx/triplet_phase_relation.h>
#include <cctbx/miller/sym_equiv.h>

namespace cctbx { namespace dmtbx {

  namespace detail {

    struct expanded_index
    {
      expanded_index(
        std::size_t ih_,
        miller::sym_equiv_index sym_equiv_index_)
      :
        ih(ih_),
        h(sym_equiv_index_.h()),
        friedel_flag(sym_equiv_index_.friedel_flag()),
        ht(sym_equiv_index_.ht())
      {}

      bool
      operator<(expanded_index const& other) const
      {
        for(std::size_t i=0;i<3;i++) {
          if (h[i] < other.h[i]) return true;
          if (h[i] > other.h[i]) return false;
        }
        return false;
      }

      std::size_t ih;
      miller::index<> h;
      bool friedel_flag;
      int ht;
    };

  } // namespace detail

  template <typename FloatType = double>
  class triplet_generator
  {
    protected:
      typedef std::map<triplet_phase_relation, std::size_t> tpr_map_type;
      typedef af::shared<tpr_map_type> list_of_tpr_maps_type;

    public:
      triplet_generator() {}

      triplet_generator(
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        bool sigma_2_only=false)
      :
        sigma_2_only_(sigma_2_only),
        t_den_(space_group.t_den())
      {
        list_of_tpr_maps_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          list_of_tpr_maps_.push_back(tpr_map_type());
        }
        std::vector<detail::expanded_index> expanded_indices;
        setup_expanded_indices(space_group, miller_indices, expanded_indices);
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          find_triplets(
            ih,
            miller_indices[ih],
            list_of_tpr_maps_[ih],
            expanded_indices);
        }
      }

      bool
      sigma_2_only() const { return sigma_2_only_; }

      int
      t_den() const { return t_den_; }

      af::shared<std::size_t>
      n_relations(bool discard_weights=false,
                  bool first_only=false) const
      {
        af::shared<std::size_t>
          result((af::reserve(list_of_tpr_maps_.size())));
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        std::size_t n_miller_indices = list_of_tpr_maps_.size();
        for(std::size_t i_h=0;i_h<n_miller_indices;i_h++,li++) {
          std::size_t prev_ik = n_miller_indices;
          std::size_t prev_ihmk = n_miller_indices;
          std::size_t n = 0;
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            if (first_only) {
              if (tpr.ik() == prev_ik && tpr.ihmk() == prev_ihmk) continue;
              prev_ik = tpr.ik();
              prev_ihmk = tpr.ihmk();
            }
            if (!discard_weights) n += lij->second;
            else n++;
          }
          result.push_back(n);
        }
        return result;
      }

      af::shared<weighted_triplet_phase_relation>
      relations_for(std::size_t ih, bool first_only=false) const
      {
        std::size_t n_miller_indices = list_of_tpr_maps_.size();
        CCTBX_ASSERT(ih < n_miller_indices);
        af::shared<weighted_triplet_phase_relation> result;
        list_of_tpr_maps_type::const_iterator
          li = list_of_tpr_maps_.begin() + ih;
        std::size_t prev_ik = n_miller_indices;
        std::size_t prev_ihmk = n_miller_indices;
        for(tpr_map_type::const_iterator
          lij=li->begin();lij!=li->end();lij++) {
          triplet_phase_relation const& tpr = lij->first;
          if (first_only) {
            if (tpr.ik() == prev_ik && tpr.ihmk() == prev_ihmk) continue;
            prev_ik = tpr.ik();
            prev_ihmk = tpr.ihmk();
          }
          result.push_back(weighted_triplet_phase_relation(tpr, lij->second));
        }
        return result;
      }

      af::shared<FloatType>
      sum_of_amplitude_products(
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<FloatType> const& amplitudes,
        bool discard_weights=false,
        bool first_only=false) const
      {
        CCTBX_ASSERT(miller_indices.size() == list_of_tpr_maps_.size());
        CCTBX_ASSERT(miller_indices.size() == amplitudes.size());
        af::shared<FloatType> result;
        result.reserve(amplitudes.size());
        list_of_tpr_maps_type::const_iterator li = list_of_tpr_maps_.begin();
        std::size_t n_miller_indices = list_of_tpr_maps_.size();
        for(std::size_t i_h=0;i_h<n_miller_indices;i_h++,li++) {
          std::size_t prev_ik = n_miller_indices;
          std::size_t prev_ihmk = n_miller_indices;
          FloatType sum = 0;
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            if (first_only) {
              if (tpr.ik() == prev_ik && tpr.ihmk() == prev_ihmk) continue;
              prev_ik = tpr.ik();
              prev_ihmk = tpr.ihmk();
            }
            FloatType a_k = amplitudes[tpr.ik()];
            FloatType a_hmk = amplitudes[tpr.ihmk()];
            FloatType aa = a_k * a_hmk;
            if (!discard_weights) aa *= lij->second;
            sum += aa;
          }
          result.push_back(sum);
        }
        return result;
      }

      af::shared<FloatType>
      apply_tangent_formula(
        af::const_ref<FloatType> const& amplitudes,
        af::const_ref<FloatType> const& phases,
        af::const_ref<bool> const& selection_fixed,
        af::const_ref<std::size_t> const& extrapolation_order,
        bool reuse_results=false,
        bool discard_weights=false,
        bool first_only=false,
        FloatType const& sum_epsilon=1.e-10) const
      {
        CCTBX_ASSERT(amplitudes.size() == list_of_tpr_maps_.size());
        CCTBX_ASSERT(phases.size() == amplitudes.size());
        CCTBX_ASSERT(   selection_fixed.size() == 0
                     || selection_fixed.size() == amplitudes.size());
        CCTBX_ASSERT(   extrapolation_order.size() == 0
                     || extrapolation_order.size() == amplitudes.size());
        CCTBX_ASSERT(first_only == false || discard_weights == true);
        af::shared<FloatType> result;
        result.assign(phases.begin(), phases.end());
        const FloatType* phase_source = (
          reuse_results ? result.begin() : phases.begin());
        std::vector<bool> fixed_or_extrapolated;
        if (selection_fixed.size() == 0) {
          fixed_or_extrapolated.resize(amplitudes.size(), false);
        }
        else {
          fixed_or_extrapolated.assign(
            selection_fixed.begin(), selection_fixed.end());
        }
        list_of_tpr_maps_type::const_iterator
          li_begin = list_of_tpr_maps_.begin();
        std::size_t i_h;
        for(std::size_t i_p=0;i_p<phases.size();i_p++) {
          if (extrapolation_order.size() == 0) {
            i_h = i_p;
          }
          else {
            i_h = extrapolation_order[i_p];
            CCTBX_ASSERT(i_h < amplitudes.size());
          }
          if (selection_fixed.size() != 0 && selection_fixed[i_h]) continue;
          CCTBX_ASSERT(!fixed_or_extrapolated[i_h]);
          list_of_tpr_maps_type::const_iterator li = li_begin + i_h;
          std::size_t prev_ik = amplitudes.size();
          std::size_t prev_ihmk = amplitudes.size();
          FloatType sum_sin(0);
          FloatType sum_cos(0);
          for(tpr_map_type::const_iterator
              lij=li->begin();lij!=li->end();lij++) {
            triplet_phase_relation const& tpr = lij->first;
            if (first_only) {
              if (tpr.ik() == prev_ik && tpr.ihmk() == prev_ihmk) continue;
              prev_ik = tpr.ik();
              prev_ihmk = tpr.ihmk();
            }
            CCTBX_ASSERT(tpr.ik() < amplitudes.size());
            CCTBX_ASSERT(tpr.ihmk() < amplitudes.size());
            if (reuse_results) {
              if (!fixed_or_extrapolated[tpr.ik()]) continue;
              if (!fixed_or_extrapolated[tpr.ihmk()]) continue;
            }
            FloatType a_k = amplitudes[tpr.ik()];
            FloatType a_hmk = amplitudes[tpr.ihmk()];
            FloatType a_k_a_hmk = a_k * a_hmk;
            if (!discard_weights) a_k_a_hmk *= lij->second;
            FloatType phi_k_phi_hmk = tpr.phi_k_phi_hmk(phase_source, t_den_);
            sum_sin += a_k_a_hmk * std::sin(phi_k_phi_hmk);
            sum_cos += a_k_a_hmk * std::cos(phi_k_phi_hmk);
          }
          if (   scitbx::fn::absolute(sum_sin) >= sum_epsilon
              || scitbx::fn::absolute(sum_cos) >= sum_epsilon) {
            result[i_h] = std::atan2(sum_sin, sum_cos);
            fixed_or_extrapolated[i_h] = true;
          }
        }
        return result;
      }

    protected:
      void
      setup_expanded_indices(
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        std::vector<detail::expanded_index>& expanded_indices)
      {
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          miller::index<> h = miller_indices[ih];
          miller::sym_equiv_indices sym_eq_h(space_group, h);
          int mult = sym_eq_h.multiplicity(false);
          for(std::size_t ih_eq=0;ih_eq<mult;ih_eq++) {
            miller::sym_equiv_index h_seq = sym_eq_h(ih_eq);
            CCTBX_ASSERT(h_seq.t_den() == t_den_);
            expanded_indices.push_back(detail::expanded_index(ih, h_seq));
          }
        }
        std::sort(expanded_indices.begin(), expanded_indices.end());
      }

      void
      find_triplets(
        std::size_t ih,
        miller::index<> const& h,
        tpr_map_type& tpr_map,
        std::vector<detail::expanded_index> const& expanded_indices)
      {
        if (expanded_indices.size() ==  0) return;
        std::size_t i_low = 0;
        std::size_t i_high = expanded_indices.size() - 1;
        const detail::expanded_index* e_low = &expanded_indices[i_low];
        const detail::expanded_index* e_high = &expanded_indices[i_high];
        while (i_low <= i_high) {
          for(std::size_t i=0;i<3;i++) {
            int s = e_low->h[i] + e_high->h[i];
            if (h[i] > s) {
              i_low++; e_low++;
              goto loop_tail;
            }
            if (h[i] < s) {
              if (i_high == 0) return;
              i_high--; e_high--;
              goto loop_tail;
            }
          }
          if (!sigma_2_only_
              || (   e_low->ih != ih
                  && e_high->ih != ih
                  && e_low->ih != e_high->ih)) {
            triplet_phase_relation tpr(
              e_low->ih,
              e_low->friedel_flag,
              e_low->ht,
              e_high->ih,
              e_high->friedel_flag,
              e_high->ht,
              t_den_);
            tpr_map[tpr] += (i_low == i_high ? 1 : 2);
          }
          i_low++; e_low++;
          i_high--; e_high--;
          loop_tail:;
        }
      }

      bool sigma_2_only_;
      int t_den_;
      list_of_tpr_maps_type list_of_tpr_maps_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_TRIPLET_GENERATOR_H
