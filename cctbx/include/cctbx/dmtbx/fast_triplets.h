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

  class miller_index_packer
  {
    public:
      miller_index_packer() {}

      miller_index_packer (af::int3 const& n) : n1_(n[1]), n2_(n[2]) {}

      long
      operator()(miller::index<> const& h) const
      {
        return (h[0] * n1_ + h[1]) * n2_ + h[2];
      }

    private:
      long n1_, n2_;
  };

  struct expanded_index_list_entry
  {
    expanded_index_list_entry(
      std::size_t ih_,
      miller::sym_equiv_index sym_equiv_index_,
      long packed_index_)
    :
      ih(ih_),
      sym_equiv_index(sym_equiv_index_),
      packed_index(packed_index_)
    {}

    bool
    operator<(expanded_index_list_entry const& other) const
    {
      return packed_index < other.packed_index;
    }

    std::size_t ih;
    miller::sym_equiv_index sym_equiv_index;
    long packed_index;
  };

  template <typename FloatType = double>
  class fast_triplets : public triplets_base<FloatType>
  {
    protected:
      typedef typename triplets_base<FloatType>::tpr_map_type tpr_map_type;
      typedef typename triplets_base<FloatType>::list_of_tpr_maps_type
                list_of_tpr_maps_type;
    private:
      typedef miller::index<> m_i;

    public:
      fast_triplets() {}

      fast_triplets(
        sgtbx::space_group_type const& sg_type,
        af::const_ref<miller::index<> > const& miller_indices)
      {
        this->t_den_ = sg_type.group().t_den();
        sgtbx::reciprocal_space::asu asu(sg_type);
        // Assert that all Miller indices are in the standard asymmetric unit.
        for(std::size_t i=0;i<miller_indices.size();i++) {
          CCTBX_ASSERT(
               miller::asym_index(sg_type.group(), asu, miller_indices[i]).h()
            == miller_indices[i]);
        }
        this->list_of_tpr_maps_.reserve(miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          this->list_of_tpr_maps_.push_back(tpr_map_type());
        }
        setup_expanded_index_list(sg_type, miller_indices);
        find_triplets(miller_indices);
      }

      void
      setup_expanded_index_list(
        sgtbx::space_group_type const& sg_type,
        af::const_ref<miller::index<> > const& miller_indices)
      {
        miller::index_span miller_index_span(miller_indices);
        //pack_index_ = miller_index_packer(miller_index_span.map_grid());
        pack_index_ = miller_index_packer(af::int3(1000,1000,1000));
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          m_i h = miller_indices[ih];
          miller::sym_equiv_indices sym_eq_h(sg_type.group(), h);
          int mult = sym_eq_h.multiplicity(false);
          for(std::size_t ih_eq=0;ih_eq<mult;ih_eq++) {
            miller::sym_equiv_index h_seq = sym_eq_h(ih_eq);
            expanded_index_list_.push_back(
              expanded_index_list_entry(ih, h_seq, pack_index_(h_seq.h())));
          }
        }
        std::sort(expanded_index_list_.begin(), expanded_index_list_.end());
        long prev_packed_index = 0;
        // sanity check
        for (std::size_t i=0;i<expanded_index_list_.size();i++) {
          long next_packed_index = expanded_index_list_[i].packed_index;
          CCTBX_ASSERT(next_packed_index != prev_packed_index);
          prev_packed_index = next_packed_index;
        }
      }

      void
      find_triplets(
        af::const_ref<miller::index<> > const& miller_indices)
      {
        for(std::size_t ih=0;ih<miller_indices.size();ih++) {
          find_triplets(ih, miller_indices[ih]);
        }
      }

      void
      find_triplets(std::size_t ih, miller::index<> const& h)
      {
        if (expanded_index_list_.size() ==  0) return;
        long unique = pack_index_(h);
        std::size_t i = 0;
        std::size_t j = expanded_index_list_.size() - 1;
        while (i <= j) {
          expanded_index_list_entry const& eile_i = expanded_index_list_[i];
          expanded_index_list_entry const& eile_j = expanded_index_list_[j];
          long sum = eile_i.packed_index + eile_j.packed_index;
          if (unique > sum) {
            i++;
          }
          else if (unique < sum) {
            if (j == 0) break;
            j--;
          }
          else {
#ifdef JUNK
            std::cout << "TPR hit"
               << " " << h.const_ref()
               << " " << eile_i.sym_equiv_index.h().const_ref()
               << " " << eile_j.sym_equiv_index.h().const_ref()
               << std::endl;
            std::cout << "       "
               << " " << unique
               << " " << eile_i.packed_index
               << " " << eile_j.packed_index
               << " " << sum
               << std::endl;
            std::cout << "       "
               << " " << (h - eile_i.sym_equiv_index.h()).const_ref()
               << std::endl;
            if (h - eile_i.sym_equiv_index.h() != eile_j.sym_equiv_index.h()) {
              std::cout << "LOOK" << std::endl;
            }
#endif
            triplet_phase_relation tpr(
              eile_i.ih,
              eile_i.sym_equiv_index,
              eile_j.ih,
              eile_j.sym_equiv_index);
            std::size_t n = (i == j ? 1 : 2);
            this->list_of_tpr_maps_[ih][tpr] += n;
            i++;
            j--;
          }
        }
      }

    private:
      std::vector<long> packed_unique_;
      miller_index_packer pack_index_;
      std::vector<expanded_index_list_entry> expanded_index_list_;
  };

}} // namespace cctbx::dmtbx

#endif // CCTBX_DMTBX_FAST_TRIPLETS_H
