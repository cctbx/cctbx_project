#ifndef CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H
#define CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H

#include <cctbx/crystal/direct_space_asu.h>

namespace cctbx { namespace crystal {

//! Next neighbor algorithms.
namespace neighbors {

  //! Simple O(N^2) algorithm for generating pairs of next neighbors.
  template <typename FloatType=double, typename IntShiftType=int>
  class simple_pair_generator
  {
    public:
      //! Convenience typedef.
      typedef typename
        direct_space_asu::asu_mappings<FloatType, IntShiftType>
          asu_mappings_t;

      //! Default constructor. Some data members are not initialized!
      simple_pair_generator() {}

      //! Initialization of the generator loop.
      /*! asu_mappings must exist and must not be changed during
          the lifetime of this instance (asu_mappings->lock()
          is called).

          distance_cutoff must be >= 0. If 0 all pairs of sites
          in the asymmetric unit and the surrounding buffer
          region are generated.
       */
      simple_pair_generator(
        asu_mappings_t* asu_mappings,
        FloatType const& distance_cutoff=0,
        bool full_matrix=false)
      :
        asu_mappings_(asu_mappings),
        distance_cutoff_sq_(distance_cutoff*distance_cutoff),
        full_matrix_(full_matrix)
      {
        CCTBX_ASSERT(distance_cutoff >= 0);
        asu_mappings->lock();
        restart();
      }

      //! Square of value as passed to the constructor.
      FloatType
      distance_cutoff_sq() const { return distance_cutoff_sq_; }

      //! Value as passed to the constructor.
      bool
      full_matrix() const { return full_matrix_; }

      /*! \brief True if the last pair returned by next() was the last
          one to be generated.
       */
      bool
      at_end() const { return at_end_; }

      //! Generates and returns the next pair.
      /*! An exception is raised if at_end() == true.
       */
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

      //! Restarts the generator.
      void
      restart()
      {
        at_end_ = false;
        incr(true);
        while (!at_end_ && pair_.dist_sq > distance_cutoff_sq_) {
          incr(false);
        }
      }

      //! Counts the number of pairs.
      std::size_t
      count_pairs()
      {
        std::size_t result = 0;
        while (!at_end_) {
          next();
          result++;
        }
        return result;
      }

      //! Calls direct_space_asu::asu_mappings::is_symmetry_interaction.
      bool
      is_symmetry_interaction(
        direct_space_asu::asu_mapping_index_pair const& pair) const
      {
        return asu_mappings_->is_symmetry_interaction(pair);
      }

    protected:
      direct_space_asu::asu_mappings<FloatType>* asu_mappings_;
      FloatType distance_cutoff_sq_;
      bool full_matrix_;
      bool at_end_;
      direct_space_asu::asu_mapping_index_pair_and_diff<FloatType> pair_;

      void
      incr(bool start);
  };

  template <typename FloatType, typename IntShiftType>
  void
  simple_pair_generator<FloatType, IntShiftType>::
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
      for(pair_.j_seq=(full_matrix_ ? 0 : pair_.i_seq);
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

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H
