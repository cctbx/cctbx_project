#ifndef CCTBX_CRYSTAL_NEIGHBORS_FAST_H
#define CCTBX_CRYSTAL_NEIGHBORS_FAST_H

#include <cctbx/crystal/direct_space_asu.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace crystal { namespace neighbors {

  template <typename FloatType=double, typename IntShiftType=int>
  class fast_pair_generator
  {
    public:
      typedef typename
        direct_space_asu::asu_mappings<FloatType, IntShiftType>
          asu_mappings_t;

      fast_pair_generator() {}

      fast_pair_generator(
        asu_mappings_t* asu_mappings,
        FloatType const& distance_cutoff,
        FloatType const& epsilon=1.e-6)
      :
        asu_mappings_(asu_mappings),
        distance_cutoff_sq_(distance_cutoff*distance_cutoff),
        epsilon_(epsilon)
      {
        CCTBX_ASSERT(distance_cutoff > 0);
        CCTBX_ASSERT(epsilon > 0);
        CCTBX_ASSERT(epsilon < 0.01);
        asu_mappings->lock();
        create_boxes(distance_cutoff * (1 + epsilon));
        restart();
      }

      FloatType
      distance_cutoff_sq() const { return distance_cutoff_sq_; }

      FloatType
      epsilon() const { return epsilon_; }

      scitbx::vec3<std::size_t> const&
      n_box() const { return n_box_; }

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
      FloatType epsilon_;
      typedef std::vector<direct_space_asu::asu_mapping_index> box_content_t;
      af::versa<box_content_t, af::c_grid<3> > boxes_;
      af::const_ref<box_content_t, af::c_grid<3> > boxes_const_ref_;
      // loop state
      bool at_end_;
      scitbx::vec3<std::size_t> n_box_;
      scitbx::vec3<std::size_t> i_box_;
      const box_content_t* boxes_i_;
      typename box_content_t::const_iterator boxes_ii_;
      scitbx::vec3<std::size_t> j_box_min_;
      scitbx::vec3<std::size_t> j_box_max_;
      scitbx::vec3<std::size_t> j_box_;
      const box_content_t* boxes_j_;
      typename box_content_t::const_iterator boxes_ji_;
      direct_space_asu::asu_mapping_index_pair_and_diff<FloatType> pair_;

      void
      create_boxes(FloatType const& box_edge)
      {
        typedef scitbx::math::float_int_conversions<FloatType, long> fic;
        cartesian<FloatType> const& box_min=asu_mappings_->mapped_sites_min();
        cartesian<FloatType> delta = asu_mappings_->mapped_sites_span();
        for(std::size_t i=0;i<3;i++) {
          n_box_[i] = static_cast<std::size_t>(
            std::max(1L, fic::iceil(delta[i] / box_edge)));
        }
        boxes_.resize(af::c_grid<3>(n_box_));
        af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
          const& mappings = asu_mappings_->mappings_const_ref();
        direct_space_asu::asu_mapping_index mi;
        for(mi.i_seq=0;mi.i_seq<mappings.size();mi.i_seq++) {
          for(mi.i_sym=0; mi.i_sym<mappings[mi.i_seq].size(); mi.i_sym++) {
            delta = mappings[mi.i_seq][mi.i_sym].mapped_site() - box_min;
            for(std::size_t i=0;i<3;i++) {
              i_box_[i] = static_cast<std::size_t>(
                fic::ifloor(delta[i] / box_edge));
              CCTBX_ASSERT(i_box_[i] < n_box_[i]);
            }
            boxes_(i_box_).push_back(mi);
          }
        }
        boxes_const_ref_ = boxes_.const_ref();
      }

      void
      incr(bool start)
      {
        af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
          const& mappings = asu_mappings_->mappings_const_ref();
        if (!start) goto continue_after_return;
        pair_.dist_sq  = -1;
        pair_.diff_vec = cartesian<FloatType>(0,0,0);
        boxes_i_ = boxes_const_ref_.begin();
        for(i_box_[0]=0;i_box_[0]<n_box_[0];i_box_[0]++) {
          j_box_min_[0] = (i_box_[0] == 0 ? 0 : i_box_[0]-1);
          j_box_max_[0] = (i_box_[0] == n_box_[0]-1 ? i_box_[0] : i_box_[0]+1);
        for(i_box_[1]=0;i_box_[1]<n_box_[1];i_box_[1]++) {
          j_box_min_[1] = (i_box_[1] == 0 ? 0 : i_box_[1]-1);
          j_box_max_[1] = (i_box_[1] == n_box_[1]-1 ? i_box_[1] : i_box_[1]+1);
        for(i_box_[2]=0;i_box_[2]<n_box_[2];i_box_[2]++, boxes_i_++) {
          j_box_min_[2] = (i_box_[2] == 0 ? 0 : i_box_[2]-1);
          j_box_max_[2] = (i_box_[2] == n_box_[2]-1 ? i_box_[2] : i_box_[2]+1);
          for(boxes_ii_=boxes_i_->begin();
              boxes_ii_!=boxes_i_->end();
              boxes_ii_++) {
            if (boxes_ii_->i_sym != 0) continue;
            pair_.i_seq = boxes_ii_->i_seq;
            for(j_box_[0]=j_box_min_[0];j_box_[0]<=j_box_max_[0];j_box_[0]++) {
            for(j_box_[1]=j_box_min_[1];j_box_[1]<=j_box_max_[1];j_box_[1]++) {
            for(j_box_[2]=j_box_min_[2];j_box_[2]<=j_box_max_[2];j_box_[2]++) {
              boxes_j_ = &boxes_const_ref_(j_box_);
              for(boxes_ji_=boxes_j_->begin();
                  boxes_ji_!=boxes_j_->end();
                  boxes_ji_++) {
                if (   boxes_ji_->i_seq <  pair_.i_seq) continue;
                if (   boxes_ji_->i_seq == pair_.i_seq
                    && boxes_ji_->i_sym == 0) continue;
                pair_.j_seq = boxes_ji_->i_seq;
                pair_.j_sym = boxes_ji_->i_sym;
                pair_.diff_vec =
                    mappings[pair_.j_seq][pair_.j_sym].mapped_site()
                  - mappings[pair_.i_seq][0].mapped_site();
                pair_.dist_sq = pair_.diff_vec.length_sq();
                return;
                continue_after_return:;
              }
            }}}
          }
        }}}
        at_end_ = true;
      }
  };

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_NEIGHBORS_FAST_H
