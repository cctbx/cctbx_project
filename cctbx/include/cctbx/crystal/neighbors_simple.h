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

      direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
      next()
      {
        CCTBX_ASSERT(!at_end_);
        direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
          result = pair_;
        incr(false);
        while (!at_end_ && pair_.dist_sq > distance_cutoff_sq_) {
          incr(false);
        }
        return result;
      }

      void
      restart()
      {
        at_end_ = false;
        incr(true);
        while (!at_end_ && pair_.dist_sq > distance_cutoff_sq_) {
          incr(false);
        }
      }

    protected:
      direct_space_asu::asu_mappings<FloatType>* asu_mappings_;
      FloatType distance_cutoff_sq_;
      direct_space_asu::asu_mapping_index_pair_and_diff<FloatType> pair_;
      std::size_t j_seq_n_sym_;
      bool at_end_;

      void
      incr(bool start)
      {
        af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
          const& mappings = asu_mappings_->mappings_const_ref();
        if (!start) goto continue_after_return;
        pair_.dist_sq  = -1;
        pair_.diff_vec = cartesian<FloatType>(0,0,0);
        for(pair_.i_seq=0;
            pair_.i_seq<mappings.size();
            pair_.i_seq++) {
          for(pair_.j_seq=pair_.i_seq;
              pair_.j_seq<mappings.size();
              pair_.j_seq++) {
            for(pair_.j_sym=(pair_.i_seq == pair_.j_seq ? 1 : 0);
                pair_.j_sym<mappings[pair_.j_seq].size();
                pair_.j_sym++) {
              if (distance_cutoff_sq_ != 0) {
                pair_.diff_vec =
                    mappings[pair_.j_seq][pair_.j_sym].mapped_site()
                  - mappings[pair_.i_seq][0].mapped_site();
                pair_.dist_sq = pair_.diff_vec.length_sq();
              }
              return;
              continue_after_return:;
            }
          }
        }
        at_end_ = true;
      }
  };

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H
