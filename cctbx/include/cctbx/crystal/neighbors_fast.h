#ifndef CCTBX_CRYSTAL_NEIGHBORS_FAST_H
#define CCTBX_CRYSTAL_NEIGHBORS_FAST_H

#include <cctbx/crystal/neighbors_simple.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace crystal { namespace neighbors {

  //! Fast algorithm for generating pairs of next neighbors.
  template <typename FloatType=double, typename IntShiftType=int>
  class fast_pair_generator
  :
    public simple_pair_generator<FloatType, IntShiftType>
  {
    public:
      //! Convenience typedef.
      typedef simple_pair_generator<FloatType, IntShiftType> base_t;
      //! Convenience typedef.
      typedef typename base_t::asu_mappings_t asu_mappings_t;

      //! Default constructor. Some data members are not initialized!
      fast_pair_generator() {}

      //! Initialization of the generator loop.
      /*! The volume of the asymmetric unit plus the buffer region
          is divided into n_boxes() boxes with edge lengths equal
          to the distance_cutoff * (1 + epsilon). The search for
          neighbors of a given site can then be restricted to the box
          of the site and the 26 neighboring boxes. If the unit cell is
          large compared to distance_cutoff this leads to a substantial
          increase in speed.

          The memory overhead for storing the boxes is roughly
          proportional to the number of sites in the asymmetric unit.

          distance_cutoff must be strictly greater than zero.

          epsilon must be greater than zero and smaller than 0.01.
       */
      fast_pair_generator(
        asu_mappings_t* asu_mappings,
        FloatType const& distance_cutoff,
        FloatType const& epsilon=1.e-6)
      :
        epsilon_(epsilon)
      {
        CCTBX_ASSERT(distance_cutoff > 0);
        CCTBX_ASSERT(epsilon > 0);
        CCTBX_ASSERT(epsilon < 0.01);
        this->asu_mappings_ = asu_mappings;
        this->distance_cutoff_sq_ = distance_cutoff*distance_cutoff;
        asu_mappings->lock();
        create_boxes(distance_cutoff * (1 + epsilon));
        restart();
      }

      //! Value as passed to the constructor.
      FloatType
      epsilon() const { return epsilon_; }

      //! Number of boxes in each dimension.
      scitbx::vec3<unsigned> const&
      n_boxes() const { return n_boxes_; }

      //! Generates and returns the next pair.
      /*! An exception is raised if at_end() == true.
       */
      direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
      next()
      {
        CCTBX_ASSERT(!this->at_end_);
        direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
          result = this->pair_;
        incr(false);
        while (  !this->at_end_
               && this->pair_.dist_sq > this->distance_cutoff_sq_) {
          incr(false);
        }
        return result;
      }

      //! Restarts the generator.
      void
      restart()
      {
        this->at_end_ = false;
        incr(true);
        while (  !this->at_end_
               && this->pair_.dist_sq > this->distance_cutoff_sq_) {
          incr(false);
        }
      }

      //! Count the number of pairs.
      std::size_t
      count_pairs()
      {
        std::size_t result = 0;
        while (!this->at_end_) {
          next();
          result++;
        }
        return result;
      }

    protected:
      FloatType epsilon_;
      typedef std::vector<direct_space_asu::asu_mapping_index> box_content_t;
      af::versa<box_content_t, af::c_grid<3, unsigned> > boxes_;
      af::const_ref<box_content_t, af::c_grid<3, unsigned> > boxes_const_ref_;
      // loop state
      scitbx::vec3<unsigned> n_boxes_;
      scitbx::vec3<unsigned> i_box_;
      const box_content_t* boxes_i_;
      typename box_content_t::const_iterator boxes_ii_;
      scitbx::vec3<unsigned> j_box_min_;
      scitbx::vec3<unsigned> j_box_max_;
      scitbx::vec3<unsigned> j_box_;
      const box_content_t* boxes_j_;
      typename box_content_t::const_iterator boxes_ji_;

      void
      create_boxes(FloatType const& box_edge);

      void
      incr(bool start);
  };

  template <typename FloatType, typename IntShiftType>
  void
  fast_pair_generator<FloatType, IntShiftType>::
  create_boxes(FloatType const& box_edge)
  {
    typedef scitbx::math::float_int_conversions<FloatType, int> fic;
    cartesian<FloatType> const&
      box_min = this->asu_mappings_->mapped_sites_min();
    cartesian<FloatType> delta = this->asu_mappings_->mapped_sites_span();
    for(std::size_t i=0;i<3;i++) {
      n_boxes_[i] = static_cast<unsigned>(
        std::max(1, fic::iceil(delta[i] / box_edge)));
    }
    boxes_.resize(af::c_grid<3, unsigned>(n_boxes_));
    af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
      const& mappings = this->asu_mappings_->mappings_const_ref();
    direct_space_asu::asu_mapping_index mi;
    for(mi.i_seq=0;mi.i_seq<mappings.size();mi.i_seq++) {
      for(mi.i_sym=0; mi.i_sym<mappings[mi.i_seq].size(); mi.i_sym++) {
        delta = mappings[mi.i_seq][mi.i_sym].mapped_site() - box_min;
        for(std::size_t i=0;i<3;i++) {
          i_box_[i] = static_cast<unsigned>(fic::ifloor(delta[i] / box_edge));
          CCTBX_ASSERT(i_box_[i] < n_boxes_[i]);
        }
        boxes_(i_box_).push_back(mi);
      }
    }
    boxes_const_ref_ = boxes_.const_ref();
  }

  template <typename FloatType, typename IntShiftType>
  void
  fast_pair_generator<FloatType, IntShiftType>::
  incr(bool start)
  {
    af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
      const& mappings = this->asu_mappings_->mappings_const_ref();
    if (!start) goto continue_after_return;
    this->pair_.dist_sq  = -1;
    this->pair_.diff_vec = cartesian<FloatType>(0,0,0);
    boxes_i_ = boxes_const_ref_.begin();
    for(i_box_[0]=0;i_box_[0]<n_boxes_[0];i_box_[0]++) {
      j_box_min_[0] = (i_box_[0] == 0 ? 0 : i_box_[0]-1);
      j_box_max_[0] = (i_box_[0] == n_boxes_[0]-1 ? i_box_[0] : i_box_[0]+1);
    for(i_box_[1]=0;i_box_[1]<n_boxes_[1];i_box_[1]++) {
      j_box_min_[1] = (i_box_[1] == 0 ? 0 : i_box_[1]-1);
      j_box_max_[1] = (i_box_[1] == n_boxes_[1]-1 ? i_box_[1] : i_box_[1]+1);
    for(i_box_[2]=0;i_box_[2]<n_boxes_[2];i_box_[2]++, boxes_i_++) {
      j_box_min_[2] = (i_box_[2] == 0 ? 0 : i_box_[2]-1);
      j_box_max_[2] = (i_box_[2] == n_boxes_[2]-1 ? i_box_[2] : i_box_[2]+1);
      for(boxes_ii_=boxes_i_->begin();
          boxes_ii_!=boxes_i_->end();
          boxes_ii_++) {
        if (boxes_ii_->i_sym != 0) continue;
        this->pair_.i_seq = boxes_ii_->i_seq;
        for(j_box_[0]=j_box_min_[0];j_box_[0]<=j_box_max_[0];j_box_[0]++) {
        for(j_box_[1]=j_box_min_[1];j_box_[1]<=j_box_max_[1];j_box_[1]++) {
        for(j_box_[2]=j_box_min_[2];j_box_[2]<=j_box_max_[2];j_box_[2]++) {
          boxes_j_ = &boxes_const_ref_(j_box_);
          for(boxes_ji_=boxes_j_->begin();
              boxes_ji_!=boxes_j_->end();
              boxes_ji_++) {
            if (   boxes_ji_->i_seq <  this->pair_.i_seq) continue;
            if (   boxes_ji_->i_seq == this->pair_.i_seq
                && boxes_ji_->i_sym == 0) continue;
            this->pair_.j_seq = boxes_ji_->i_seq;
            this->pair_.j_sym = boxes_ji_->i_sym;
            this->pair_.diff_vec =
               mappings[this->pair_.j_seq][this->pair_.j_sym].mapped_site()
             - mappings[this->pair_.i_seq][0].mapped_site();
            this->pair_.dist_sq = this->pair_.diff_vec.length_sq();
            return;
            continue_after_return:;
          }
        }}}
      }
    }}}
    this->at_end_ = true;
  }

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_NEIGHBORS_FAST_H
