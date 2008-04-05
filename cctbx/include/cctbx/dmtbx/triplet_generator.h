/* References:
     C.M. Weeks, P.D. Adams, J. Berendzen, A.T. Brunger, E.J. Dodson,
     R.W. Grosse-Kunstleve, T.R. Schneider, G.M. Sheldrick,
     T.C. Terwilliger, M. Turkenburg, I. Uson
     Automatic solution of heavy-atom substructures.
     Methods in Enzymology 2003, 374, 37-83.

     Sheldrick, G.M. (1982).
     Crystallographic algorithms for mini- and maxi-computers.
     In: Computational Crystallography, Ed. D. Sayre,
     Oxford University Press, 506-514.
 */

#ifndef CCTBX_DMTBX_TRIPLET_GENERATOR_H
#define CCTBX_DMTBX_TRIPLET_GENERATOR_H

#include <cctbx/dmtbx/triplet_phase_relation.h>
#include <cctbx/miller/sym_equiv.h>
#include <scitbx/array_family/sort.h>
#include <scitbx/array_family/selections.h>
#include <boost/scoped_array.hpp>

namespace cctbx {

//! Direct methods toolbox.
namespace dmtbx {

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

  //! Triplets for direct methods (reciprocal-space squaring).
  template <typename FloatType = double>
  class triplet_generator
  {
    protected:
      typedef weighted_triplet_phase_relation wtpr_t;
      typedef af::shared<af::shared<wtpr_t> > array_of_wtprs_t;
      typedef af::const_ref<wtpr_t> cr_wtprs_t;

    public:
      //! Default constructor. Some data members are not initialized!
      triplet_generator() {}

      //! Searches for all triplets h = k + h-k.
      /*! A triplet phase relation h = k + h-k is a "sigma-2" relation
          only if the three Miller indices involved are not related
          by symmetry.

          If amplitudes are given and max_relations_per_reflection > 0
          only the triplets with the max_relations_per_reflection largest
          values of the product amplitudes[k]*amplitudes[h-k]
          are retained.

          By default this class keeps track of the number of times
          each triplet is found. This can be disabled with
          discard_weights=true.

          With sigmas_2_only=false and discard_weights=false application
          of the tangent formula is exactly equivalent to squaring
          in direct space.
       */
      triplet_generator(
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<FloatType> const&
          amplitudes=af::const_ref<FloatType>(0,0),
        std::size_t max_relations_per_reflection=0,
        bool sigma_2_only=false,
        bool discard_weights=false)
      :
        t_den_(space_group.t_den()),
        max_relations_per_reflection_(max_relations_per_reflection),
        sigma_2_only_(sigma_2_only),
        discard_weights_(discard_weights),
        array_of_wtprs_((af::reserve(miller_indices.size())))
      {
        CCTBX_ASSERT(   amplitudes.size() == 0
                     || amplitudes.size() == miller_indices.size());
        CCTBX_ASSERT(   max_relations_per_reflection == 0
                     || amplitudes.size() > 0);
        std::vector<detail::expanded_index> expanded_indices;
        setup_expanded_indices(space_group, miller_indices, expanded_indices);
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          af::shared<wtpr_t> tprs = find_triplets(
            ih,
            miller_indices[ih],
            expanded_indices);
          if (   max_relations_per_reflection == 0
              || max_relations_per_reflection >= tprs.size()) {
            array_of_wtprs_.push_back(tprs);
          }
          else {
            array_of_wtprs_.push_back(truncate(
              tprs.const_ref(),
              amplitudes,
              max_relations_per_reflection));
          }
        }
      }

      //! Translation part denominator of the space_group constructor argument.
      int
      t_den() const { return t_den_; }

      //! Value of max_relations_per_reflection_ passed to the constructor.
      std::size_t
      max_relations_per_reflection() const
      {
        return max_relations_per_reflection_;
      }

      //! Value of sigma_2_only flag passed to the constructor.
      bool
      sigma_2_only() const { return sigma_2_only_; }

      //! Value of discard_weights flag passed to the constructor.
      bool
      discard_weights() const { return discard_weights_; }

      //! Number of triplet phase relations for each unique Miller index.
      /*! The Miller indices are implied by the position in the array.
       */
      af::shared<std::size_t>
      n_relations() const
      {
        af::shared<std::size_t> result((af::reserve(array_of_wtprs_.size())));
        std::size_t n_miller_indices = array_of_wtprs_.size();
        for(std::size_t ih=0;ih<n_miller_indices;ih++) {
          cr_wtprs_t tprs = array_of_wtprs_[ih].const_ref();
          std::size_t n = 0;
          for(const wtpr_t* tpr=tprs.begin();tpr!=tprs.end();tpr++) {
            n += tpr->weight();
          }
          result.push_back(n);
        }
        return result;
      }

      //! Array of triplet phase relations for a selected Miller index.
      /*! The Miller indices are implied by the position in the array.
       */
      af::shared<weighted_triplet_phase_relation>
      relations_for(std::size_t ih)
      {
        std::size_t n_miller_indices = array_of_wtprs_.size();
        CCTBX_ASSERT(ih < n_miller_indices);
        return array_of_wtprs_[ih];
      }

      //! Sum of amplitude products for each unique Miller index.
      /*! The Miller indices are implied by the position in the array.
       */
      af::shared<FloatType>
      sums_of_amplitude_products(
        af::const_ref<FloatType> const& amplitudes) const
      {
        std::size_t n_miller_indices = array_of_wtprs_.size();
        CCTBX_ASSERT(amplitudes.size() == n_miller_indices);
        af::shared<FloatType> result((af::reserve(n_miller_indices)));
        for(std::size_t ih=0;ih<n_miller_indices;ih++) {
          cr_wtprs_t tprs = array_of_wtprs_[ih].const_ref();
          FloatType sum = 0;
          for(const wtpr_t* tpr=tprs.begin();tpr!=tprs.end();tpr++) {
            sum += amplitudes[tpr->ik()]
                 * amplitudes[tpr->ihmk()]
                 * tpr->weight();
          }
          result.push_back(sum);
        }
        return result;
      }

      //! Computation of new phase estimates.
      /*! The Miller indices are implied by the position in the array.

          Phases for which selection_fixed is true are not changed.

          If use_fixed_only is true only the fixed phases are
          used in the computation of the new phases. Otherwise
          all phases are used.

          If reuse_results is true newly computed phases are
          used subsequently in the computation of other phases.

          If both use_fixed_only and reuse_results are true,
          both fixed phases and newly computed phases are
          used in the computation of the new phases.
       */
      af::shared<FloatType>
      apply_tangent_formula(
        af::const_ref<FloatType> const& amplitudes,
        af::const_ref<FloatType> const& phases_rad,
        af::const_ref<bool> const& selection_fixed,
        bool use_fixed_only=false,
        bool reuse_results=false,
        FloatType const& sum_epsilon=1.e-10) const
      {
        CCTBX_ASSERT(amplitudes.size() == array_of_wtprs_.size());
        CCTBX_ASSERT(phases_rad.size() == amplitudes.size());
        CCTBX_ASSERT(   selection_fixed.size() == 0
                     || selection_fixed.size() == amplitudes.size());
        CCTBX_ASSERT(!use_fixed_only || selection_fixed.size() > 0);
        af::shared<FloatType> result(phases_rad.begin(), phases_rad.end());
        const FloatType* phase_source = (
          reuse_results ? result.begin() : phases_rad.begin());
        std::vector<bool> fixed_or_extrapolated;
        if (selection_fixed.size() == 0) {
          fixed_or_extrapolated.resize(amplitudes.size(), false);
        }
        else {
          fixed_or_extrapolated.assign(
            selection_fixed.begin(), selection_fixed.end());
        }
        for(std::size_t ih=0;ih<phases_rad.size();ih++) {
          if (selection_fixed.size() != 0 && selection_fixed[ih]) continue;
          CCTBX_ASSERT(!fixed_or_extrapolated[ih]);
          cr_wtprs_t tprs = array_of_wtprs_[ih].const_ref();
          FloatType sum_sin(0);
          FloatType sum_cos(0);
          for(const wtpr_t* tpr=tprs.begin();tpr!=tprs.end();tpr++) {
            CCTBX_ASSERT(tpr->ik() < amplitudes.size());
            CCTBX_ASSERT(tpr->ihmk() < amplitudes.size());
            if (use_fixed_only) {
              if (reuse_results) {
                if (!fixed_or_extrapolated[tpr->ik()]) continue;
                if (!fixed_or_extrapolated[tpr->ihmk()]) continue;
              }
              else {
                if (!selection_fixed[tpr->ik()]) continue;
                if (!selection_fixed[tpr->ihmk()]) continue;
              }
            }
            FloatType a_k_a_hmk = amplitudes[tpr->ik()]
                                * amplitudes[tpr->ihmk()]
                                * tpr->weight();
            FloatType phi_k_phi_hmk = tpr->phi_k_phi_hmk(phase_source, t_den_);
            sum_sin += a_k_a_hmk * std::sin(phi_k_phi_hmk);
            sum_cos += a_k_a_hmk * std::cos(phi_k_phi_hmk);
          }
          if (   scitbx::fn::absolute(sum_sin) >= sum_epsilon
              || scitbx::fn::absolute(sum_cos) >= sum_epsilon) {
            result[ih] = std::atan2(sum_sin, sum_cos);
            fixed_or_extrapolated[ih] = true;
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

      struct expanded_indices_scanner
      {
        expanded_indices_scanner(
          std::vector<detail::expanded_index> const& expanded_indices)
        :
          i_low(0),
          i_high(expanded_indices.size() - 1),
          e_low(&expanded_indices[i_low]),
          e_high(&expanded_indices[i_high])
        {}

        bool
        incr_low()
        {
          if (i_low == i_high) return false;
          i_low++;
          e_low++;
          return true;
        }

        bool
        decr_high()
        {
          if (i_low == i_high) return false;
          i_high--;
          e_high--;
          return true;
        }

        bool
        advance()
        {
          if (!incr_low()) return false;
          if (!decr_high()) return false;
          return true;
        }

        bool
        find_next(miller::index<> const& h)
        {
          for(std::size_t i=0;i<3;) {
            int s = e_low->h[i] + e_high->h[i];
            if (h[i] > s) {
              if (!incr_low()) return false;
              i = 0;
            }
            else if (h[i] < s) {
              if (!decr_high()) return false;
              i = 0;
            }
            else {
              i++;
            }
          }
          return true;
        }

        bool
        current_is_sigma_2(std::size_t ih) const
        {
          return e_low->ih != ih
              && e_high->ih != ih
              && e_low->ih != e_high->ih;
        }

        triplet_phase_relation
        get_tpr(int t_den) const
        {
          return triplet_phase_relation(
            e_low->ih,
            e_low->friedel_flag,
            e_low->ht,
            e_high->ih,
            e_high->friedel_flag,
            e_high->ht,
            t_den);
        }

        std::size_t
        get_weight() const
        {
          if (i_low == i_high) return 1;
          return 2;
        }

        std::size_t i_low;
        std::size_t i_high;
        const detail::expanded_index* e_low;
        const detail::expanded_index* e_high;
      };

      af::shared<weighted_triplet_phase_relation>
      find_triplets(
        std::size_t ih,
        miller::index<> const& h,
        std::vector<detail::expanded_index> const& expanded_indices)
      {
        typedef std::map<triplet_phase_relation, std::size_t> tpr_map_t;
        tpr_map_t tpr_map;
        tpr_map_t::const_iterator m;
        if (expanded_indices.size() != 0) {
          expanded_indices_scanner scanner(expanded_indices);
          while (scanner.find_next(h)) {
            if (!sigma_2_only_ || scanner.current_is_sigma_2(ih)) {
              tpr_map[scanner.get_tpr(t_den_)] += scanner.get_weight();
            }
            if (!scanner.advance()) break;
          }
        }
        af::shared<wtpr_t> wtpr_array((af::reserve(tpr_map.size())));
        if (!discard_weights_) {
          for(m=tpr_map.begin();m!=tpr_map.end();m++) {
            wtpr_array.push_back(wtpr_t(m->first, m->second));
          }
        }
        else {
          const triplet_phase_relation* prev_tpr = 0;
          for(m=tpr_map.begin();m!=tpr_map.end();m++) {
            if (prev_tpr != 0 && m->first.is_similar_to(*prev_tpr)) continue;
            prev_tpr = &m->first;
            wtpr_array.push_back(wtpr_t(m->first, 1));
          }
        }
        return wtpr_array;
      }

      af::shared<weighted_triplet_phase_relation>
      truncate(
        af::const_ref<weighted_triplet_phase_relation> const& tprs,
        af::const_ref<FloatType> const& amplitudes,
        std::size_t max_relations_per_reflection)
      {
        CCTBX_ASSERT(tprs.size() > max_relations_per_reflection);
        af::shared<weighted_triplet_phase_relation>
          result((af::reserve(max_relations_per_reflection)));
        boost::scoped_array<FloatType>
          values_to_sort(new FloatType[tprs.size()]);
        FloatType* v = values_to_sort.get();
        for(const wtpr_t* tpr=tprs.begin();tpr!=tprs.end();tpr++) {
          *v++ = amplitudes[tpr->ik()]
               * amplitudes[tpr->ihmk()]
               * tpr->weight();
        }
        af::shared<std::size_t> perm = af::sort_permutation(
          af::const_ref<FloatType>(values_to_sort.get(), tprs.size()),
          true);
        perm.resize(max_relations_per_reflection);
        return af::select(tprs, perm.const_ref());
      }

      int t_den_;
      std::size_t max_relations_per_reflection_;
      bool sigma_2_only_;
      bool discard_weights_;
      array_of_wtprs_t array_of_wtprs_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_TRIPLET_GENERATOR_H
