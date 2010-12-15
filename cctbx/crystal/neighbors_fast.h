#ifndef CCTBX_CRYSTAL_NEIGHBORS_FAST_H
#define CCTBX_CRYSTAL_NEIGHBORS_FAST_H

#include <cctbx/crystal/neighbors_simple.h>
#include <scitbx/cubicles.h>
#include <map>
#include <set>

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

          If minimal == true pairs with i_seq > j_seq will be
          suppressed even if j_sym != 0.
          See also:
            cctbx::crystal::direct_space_asu::asu_mapping_index_pair::is_active

          epsilon must be greater than zero and smaller than 0.01.
       */
      fast_pair_generator(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<
            FloatType, IntShiftType> > const& asu_mappings,
        FloatType const& distance_cutoff,
        bool minimal=false,
        FloatType const& min_cubicle_edge=5,
        FloatType const& epsilon=1e-6)
      :
        epsilon_(epsilon),
        cubicles_(
          asu_mappings.get()->mapped_sites_min(),
          asu_mappings.get()->mapped_sites_span(),
          std::max(distance_cutoff, min_cubicle_edge),
          epsilon),
        n_boxes_(cubicles_.ref.accessor())
      {
        CCTBX_ASSERT(epsilon > 0);
        CCTBX_ASSERT(epsilon < 0.01);
        this->asu_mappings_owner_ = asu_mappings;
        this->asu_mappings_ = asu_mappings.get();
        this->distance_cutoff_sq_ = distance_cutoff*distance_cutoff;
        this->minimal_ = minimal;
        af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
          const& mappings = this->asu_mappings_->mappings_const_ref();
        direct_space_asu::asu_mapping_index mi;
        for(mi.i_seq=0;mi.i_seq<mappings.size();mi.i_seq++) {
          for(mi.i_sym=0; mi.i_sym<mappings[mi.i_seq].size(); mi.i_sym++) {
            std::size_t i1d_cub = cubicles_.ref.accessor()(
              cubicles_.i_cubicle(mappings[mi.i_seq][mi.i_sym].mapped_site()));
            cubicles_.ref[i1d_cub].push_back(mi);
          }
        }
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

      //! Counts the number of pairs.
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
                  == this->asu_mappings_->mappings_const_ref().size());
        af::shared<bool> result(
          primary_selection.begin(),
          primary_selection.end());
        af::ref<bool> result_ = result.ref();
        while (!this->at_end_) {
          direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
            pair = next();
          if      (primary_selection[pair.i_seq]) result_[pair.j_seq] = true;
          else if (primary_selection[pair.j_seq]) result_[pair.i_seq] = true;
        }
        return result;
      }

      af::shared< std::set< unsigned > >
      distance_based_simple_two_way_bond_sets(
        af::const_ref< std::string > const& elements,
        af::const_ref<std::size_t> const& conformer_indices,
        std::map<std::string, double> expected_bond_lengths,
        std::map<std::string, double> vdw_radii,
        double fallback_expected_bond_length,
        double tolerance_factor_expected_bond_length)
      {
        CCTBX_ASSERT(conformer_indices.size() == elements.size());
        std::size_t n_sites = elements.size();
        af::shared< std::set< unsigned > > bonds(n_sites);
        restart();
        while (!this->at_end_) {
          direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
            pair = next();
          std::size_t conf1 = conformer_indices[pair.i_seq];
          std::size_t conf2 = conformer_indices[pair.j_seq];
          if ((conf2 != conf1) && (conf1 != 0) && (conf2 != 0)) {
            continue;
          }
          std::string const& elem1 = elements[pair.i_seq];
          std::string const& elem2 = elements[pair.j_seq];
          std::string elem_key;
          if (elem1 < elem2) {
            elem_key = elem1 + ":" + elem2;
          } else {
            elem_key = elem2 + ":" + elem1;
          }
          std::map<std::string, double>::const_iterator
            ebl_from_tab = expected_bond_lengths.find(elem_key);
          double ebl;
          if (ebl_from_tab != expected_bond_lengths.end()) {
            ebl = ebl_from_tab->second;
            if (ebl == 0.0) continue;
          } else {
            double radius1 = 0.0;
            double radius2 = 0.0;
            if (vdw_radii.find(elem1) != vdw_radii.end()) {
              radius1 = vdw_radii[elem1];
            }
            if (vdw_radii.find(elem2) != vdw_radii.end()) {
              radius2 = vdw_radii[elem2];
            }
            if (radius1 > radius2) {
              ebl = radius1;
            } else {
              ebl = radius2;
            }
            if (ebl == 0.0) {
              ebl = fallback_expected_bond_length;
            }
            if (ebl < 0.0) continue;
          }
          double cutoff = ebl * tolerance_factor_expected_bond_length;
          double cutoff_sq = cutoff * cutoff;
          if (pair.dist_sq > cutoff_sq) continue;
          bonds[pair.i_seq].insert(pair.j_seq);
          bonds[pair.j_seq].insert(pair.i_seq);
        }
        return bonds;
      }


    protected:
      FloatType epsilon_;
      typedef std::vector<direct_space_asu::asu_mapping_index> box_content_t;
      scitbx::cubicles<box_content_t, FloatType> cubicles_;
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
      incr(bool start);

      bool
      is_active_pair(
        unsigned i_seq,
        direct_space_asu::asu_mapping_index const& other) const
      {
        if (other.i_seq >  i_seq) return true;
        if (other.i_seq == i_seq) return (other.i_sym != 0);
        return (!this->minimal_ && other.i_sym != 0);
      }
  };

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
    boxes_i_ = cubicles_.ref.begin();
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
        for(j_box_[0]=j_box_min_[0];j_box_[0]<=j_box_max_[0];j_box_[0]++) {
        for(j_box_[1]=j_box_min_[1];j_box_[1]<=j_box_max_[1];j_box_[1]++) {
        for(j_box_[2]=j_box_min_[2];j_box_[2]<=j_box_max_[2];j_box_[2]++) {
          boxes_j_ = &cubicles_.ref(j_box_);
          for(boxes_ji_=boxes_j_->begin();
              boxes_ji_!=boxes_j_->end();
              boxes_ji_++) {
            if (!is_active_pair(boxes_ii_->i_seq, *boxes_ji_)) continue;
            this->pair_.i_seq = boxes_ii_->i_seq;
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
