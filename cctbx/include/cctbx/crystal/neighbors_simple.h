#ifndef CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H
#define CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace crystal { namespace neighbors {

  template <typename FloatType=double, typename IntShiftType=int>
  class simple_pair_generator
  {
    public:
      typedef typename
        direct_space_asu::asu_mappings<FloatType, IntShiftType>
          asu_mappings_t;

      simple_pair_generator() {}

      simple_pair_generator(
        asu_mappings_t* asu_mappings,
        FloatType const& distance_cutoff=0)
      :
        asu_mappings_(asu_mappings),
        distance_cutoff_sq_(distance_cutoff*distance_cutoff)
      {
        asu_mappings->lock();
        restart();
      }

      bool
      at_end() const { return at_end_; }

      direct_space_asu::asu_mapping_index_pair<FloatType>
      next()
      {
        CCTBX_ASSERT(!at_end_);
        direct_space_asu::asu_mapping_index_pair<FloatType> result = pair_;
        incr();
        while (!at_end_ && pair_.dist_sq > distance_cutoff_sq_) {
          incr();
        }
        return result;
      }

      void
      restart()
      {
        at_end_ = false;
        incr_to_first();
        while (!at_end_ && pair_.dist_sq > distance_cutoff_sq_) {
          incr();
        }
      }

    protected:
      direct_space_asu::asu_mappings<FloatType>* asu_mappings_;
      FloatType distance_cutoff_sq_;
      direct_space_asu::asu_mapping_index_pair<FloatType> pair_;
      std::size_t j_seq_n_sym_;
      bool at_end_;

      typedef
        af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
          mappings_ref_t;

      void
      incr_to_first()
      {
        pair_.dist_sq  = -1;
        mappings_ref_t mappings = asu_mappings_->mappings().const_ref();
        for(pair_.i_seq=0;
            pair_.i_seq<mappings.size();
            pair_.i_seq++) {
          for(pair_.j_seq=pair_.i_seq;
              pair_.j_seq<mappings.size();
              pair_.j_seq++) {
            pair_.j_sym = 0;
            j_seq_n_sym_ = mappings[pair_.j_seq].size();
            if (pair_.i_seq == pair_.j_seq) {
              if (j_seq_n_sym_ > 1) {
                pair_.j_sym++;
                set_pair_dist_sq(mappings);
                return;
              }
            }
            else if (j_seq_n_sym_ > 0) {
              set_pair_dist_sq(mappings);
              return;
            }
          }
        }
        at_end_ = true;
      }

      void
      incr()
      {
        mappings_ref_t mappings = asu_mappings_->mappings().const_ref();
        pair_.j_sym++;
        while (pair_.j_sym == j_seq_n_sym_) {
          pair_.j_sym = 0;
          pair_.j_seq++;
          while (pair_.j_seq == mappings.size()) {
            pair_.i_seq++;
            if (pair_.i_seq == mappings.size()) {
              pair_.i_seq = 0;
              pair_.dist_sq = -1;
              at_end_ = true;
              return;
            }
            pair_.j_seq = pair_.i_seq;
            j_seq_n_sym_ = mappings[pair_.j_seq].size();
            if (j_seq_n_sym_ > 1) {
              pair_.j_sym++;
              set_pair_dist_sq(mappings);
              return;
            }
            pair_.j_seq++;
          }
          j_seq_n_sym_ = mappings[pair_.j_seq].size();
        }
        set_pair_dist_sq(mappings);
      }

      void
      set_pair_dist_sq(mappings_ref_t const& mappings)
      {
        if (distance_cutoff_sq_ == 0) return;
        pair_.dist_sq = (  mappings[pair_.j_seq][pair_.j_sym].mapped_site()
                         - mappings[pair_.i_seq][0].mapped_site()).length_sq();
      }
  };

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H
