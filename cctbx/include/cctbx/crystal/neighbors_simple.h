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
      /*! distance_cutoff must be >= 0. If 0 all pairs of sites
          in the asymmetric unit and the surrounding buffer
          region are generated.

          If minimal == true pairs with i_seq > j_seq will be
          suppressed even if j_sym != 0.
          See also:
            cctbx::crystal::direct_space_asu::asu_mapping_index_pair::is_active
       */
      simple_pair_generator(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<
            FloatType, IntShiftType> > const& asu_mappings,
        FloatType const& distance_cutoff=0,
        bool minimal=false)
      :
        asu_mappings_owner_(asu_mappings),
        asu_mappings_(asu_mappings.get()),
        distance_cutoff_sq_(distance_cutoff*distance_cutoff),
        minimal_(minimal)
      {
        CCTBX_ASSERT(distance_cutoff >= 0);
        restart();
      }

      //! Instance as passed to the constructor.
      boost::shared_ptr<
        direct_space_asu::asu_mappings<
          FloatType, IntShiftType> >
      asu_mappings() const { return asu_mappings_owner_; }

      //! Square of value as passed to the constructor.
      FloatType
      distance_cutoff_sq() const { return distance_cutoff_sq_; }

      //! Value as passed to the constructor.
      bool
      minimal() const { return minimal_; }

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
        is_swapped_ = false;
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

      //! Maximum distance squared of all remaining pairs.
      FloatType
      max_distance_sq()
      {
        FloatType result = -1;
        while (!this->at_end_) {
          result = std::max(result, next().dist_sq);
        }
        return result;
      }

      /*! \brief Selection of all neighbors within distance_cutoff
          of primary_selection.
       */
      /*! The result includes the primary_selection.
       */
      af::shared<bool>
      neighbors_of(af::const_ref<bool> const& primary_selection)
      {
        CCTBX_ASSERT(primary_selection.size()
                  == asu_mappings_->mappings_const_ref().size());
        af::shared<bool> result(
          primary_selection.begin(),
          primary_selection.end());
        af::ref<bool> result_ = result.ref();
        while (!at_end_) {
          direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
            pair = next();
          if      (primary_selection[pair.i_seq]) result_[pair.j_seq] = true;
          else if (primary_selection[pair.j_seq]) result_[pair.i_seq] = true;
        }
        return result;
      }

      //! Calls direct_space_asu::asu_mappings::is_simple_interaction.
      bool
      is_simple_interaction(
        direct_space_asu::asu_mapping_index_pair const& pair) const
      {
        return asu_mappings_->is_simple_interaction(pair);
      }

    protected:
      boost::shared_ptr<asu_mappings_t> asu_mappings_owner_;
      const asu_mappings_t* asu_mappings_;
      unsigned mappings_size_;
      FloatType distance_cutoff_sq_;
      bool minimal_;
      bool at_end_;
      bool is_swapped_;
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
    if (!start) {
      if (!is_swapped_) goto continue_after_return;
      goto continue_after_return_swapped;
    }
    mappings_size_ = static_cast<unsigned>(mappings.size());
    pair_.dist_sq  = -1;
    pair_.diff_vec = cartesian<FloatType>(0,0,0);
    for(pair_.i_seq=0;
        pair_.i_seq<mappings_size_;
        pair_.i_seq++) {
      for(pair_.j_seq=pair_.i_seq;
          pair_.j_seq<mappings_size_;
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
        if (!minimal_ && pair_.j_seq != pair_.i_seq) {
          std::swap(pair_.i_seq, pair_.j_seq);
          is_swapped_ = true;
          for(pair_.j_sym=1;
              pair_.j_sym<mappings[pair_.j_seq].size();
              pair_.j_sym++) {
            if (distance_cutoff_sq_ != 0) {
              pair_.diff_vec =
                  mappings[pair_.j_seq][pair_.j_sym].mapped_site()
                - mappings[pair_.i_seq][0].mapped_site();
              pair_.dist_sq = pair_.diff_vec.length_sq();
            }
            return;
            continue_after_return_swapped:;
          }
          std::swap(pair_.i_seq, pair_.j_seq);
          is_swapped_ = false;
        }
      }
    }
    at_end_ = true;
  }

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_NEIGHBORS_SIMPLE_H
