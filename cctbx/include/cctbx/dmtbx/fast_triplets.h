/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created (R.W. Grosse-Kunstleve)
 */

/* Reference:
     Sheldrick, G.M. (1982).
     Crystallographic algorithms for mini- and maxi-computers.
     In: Computational Crystallography, Ed. D. Sayre,
     Oxford University Press, 506-514.
 */

#ifndef CCTBX_DMTBX_FAST_TRIPLETS_H
#define CCTBX_DMTBX_FAST_TRIPLETS_H

#include <cctbx/dmtbx/triplet.h>

namespace cctbx { namespace dmtbx {

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

  template <typename FloatType = double>
  class fast_triplets : public triplets_base<FloatType>
  {
    protected:
      typedef typename triplets_base<FloatType>::tpr_map_type tpr_map_type;
      typedef typename triplets_base<FloatType>::list_of_tpr_maps_type
                list_of_tpr_maps_type;

    public:
      fast_triplets() {}

      fast_triplets(
        sgtbx::space_group_type const& sg_type,
        af::const_ref<miller::index<> > const& miller_indices,
        bool sigma_2_only=false)
      :
        sigma_2_only_(sigma_2_only)
      {
        this->t_den_ = sg_type.group().t_den();
        this->list_of_tpr_maps_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          this->list_of_tpr_maps_.push_back(tpr_map_type());
        }
        std::vector<expanded_index> expanded_indices;
        setup_expanded_indices(sg_type, miller_indices, expanded_indices);
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          find_triplets(
            ih,
            miller_indices[ih],
            this->list_of_tpr_maps_[ih],
            expanded_indices);
        }
      }

      bool
      sigma_2_only() const { return sigma_2_only_; }

    protected:
      void
      setup_expanded_indices(
        sgtbx::space_group_type const& sg_type,
        af::const_ref<miller::index<> > const& miller_indices,
        std::vector<expanded_index>& expanded_indices)
      {
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          miller::index<> h = miller_indices[ih];
          miller::sym_equiv_indices sym_eq_h(sg_type.group(), h);
          int mult = sym_eq_h.multiplicity(false);
          for(std::size_t ih_eq=0;ih_eq<mult;ih_eq++) {
            miller::sym_equiv_index h_seq = sym_eq_h(ih_eq);
            CCTBX_ASSERT(h_seq.t_den() == this->t_den_);
            expanded_indices.push_back(expanded_index(ih, h_seq));
          }
        }
        std::sort(expanded_indices.begin(), expanded_indices.end());
      }

      void
      find_triplets(
        std::size_t ih,
        miller::index<> const& h,
        tpr_map_type& tpr_map,
        std::vector<expanded_index> const& expanded_indices)
      {
        if (expanded_indices.size() ==  0) return;
        std::size_t i_low = 0;
        std::size_t i_high = expanded_indices.size() - 1;
        const expanded_index* e_low = &expanded_indices[i_low];
        const expanded_index* e_high = &expanded_indices[i_high];
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
              this->t_den_);
            tpr_map[tpr] += (i_low == i_high ? 1 : 2);
          }
          i_low++; e_low++;
          i_high--; e_high--;
          loop_tail:;
        }
      }

    private:
      bool sigma_2_only_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_FAST_TRIPLETS_H
