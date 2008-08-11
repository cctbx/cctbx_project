#ifndef CCTBX_CRYSTAL_ASU_CLUSTERS_H
#define CCTBX_CRYSTAL_ASU_CLUSTERS_H

#include <cctbx/crystal/pair_tables.h>
#include <scitbx/array_family/sort.h>

namespace cctbx { namespace crystal {

  /*! \brief Determination of clusters of points with interactions
      defined by pair_asu_table.
   */
  class asu_clusters
  {
    public:
      //! Default constructor.
      asu_clusters() {}

      //! Initialization of index_groups.
      /*! The index_groups are not sorted by size.
          The indices in each group are also not sorted.
       */
      template <typename FloatType, typename IntShiftType>
      asu_clusters(
        crystal::pair_asu_table<FloatType, IntShiftType> const& pair_asu_table,
        bool strictly_in_asu=true)
      {
        af::const_ref<pair_asu_dict> table=pair_asu_table.table().const_ref();
        unsigned n_sites = table.size();
        std::vector<unsigned> i_clusters(n_sites, n_sites);
        std::vector<std::vector<unsigned> > index_groups_;
        unsigned n_cleared_index_groups = 0;
        for(unsigned i_seq=0;i_seq<n_sites;i_seq++) {
          unsigned i_cluster = i_clusters[i_seq];
          if (i_cluster == n_sites) {
            i_cluster = i_clusters[i_seq] = index_groups_.size();
            index_groups_.push_back(std::vector<unsigned>(1, i_seq));
          }
          std::vector<unsigned>* index_group_i = &index_groups_[i_cluster];
          for(pair_asu_dict::const_iterator
                asu_dict_i = table[i_seq].begin();
                asu_dict_i != table[i_seq].end();
                asu_dict_i++) {
            unsigned j_seq = asu_dict_i->first;
            unsigned j_cluster = i_clusters[j_seq];
            if (j_cluster == i_cluster) continue;
            for(pair_asu_j_sym_groups::const_iterator
                  j_sym_groups_i = asu_dict_i->second.begin();
                  j_sym_groups_i != asu_dict_i->second.end();
                  j_sym_groups_i++) {
              unsigned j_syms_0 = *(j_sym_groups_i->begin());
              if (!strictly_in_asu || j_syms_0 == 0) {
                if (j_cluster == n_sites) {
                  i_clusters[j_seq] = i_cluster;
                  index_group_i->push_back(j_seq);
                }
                else {
                  for(unsigned i=0;i<index_group_i->size();i++) {
                    i_clusters[(*index_group_i)[i]] = j_cluster;
                  }
                  std::vector<unsigned>*
                    index_groups_j = &index_groups_[j_cluster];
                  index_groups_j->reserve(
                    index_groups_j->size() + index_group_i->size());
                  for(unsigned i=0;i<index_group_i->size();i++) {
                    index_groups_j->push_back((*index_group_i)[i]);
                  }
                  index_group_i = index_groups_j;
                  // free memory
                  index_groups_[i_cluster] = std::vector<unsigned>();
                  n_cleared_index_groups++;
                  i_cluster = j_cluster;
                }
                break;
              }
            }
          }
        }
        index_groups.resize(index_groups_.size() - n_cleared_index_groups);
        af::ref<std::vector<unsigned> > index_groups_ref = index_groups.ref();
        std::size_t j_cluster=0;
        for(std::size_t
              i_cluster=0;
              i_cluster<index_groups_.size();
              i_cluster++) {
          if (index_groups_[i_cluster].size() > 0) {
            // swap is used to avoid deep copies of index_groups
            index_groups_ref[j_cluster++].swap(index_groups_[i_cluster]);
          }
        }
        CCTBX_ASSERT(j_cluster == index_groups_ref.size());
      }

      //! In-place sorting of index_groups.
      asu_clusters&
      sort_index_groups_by_size()
      {
        af::ref<std::vector<unsigned> > index_groups_ref = index_groups.ref();
        std::vector<unsigned> sizes;
        sizes.reserve(index_groups_ref.size());
        for(std::size_t i=0;i<index_groups_ref.size();i++) {
          sizes.push_back(index_groups_ref[i].size());
        }
        af::shared<std::size_t> perm_owner = af::sort_permutation(
          af::make_ref(sizes), true);
        af::const_ref<std::size_t> perm = perm_owner.const_ref();
        af::shared<std::vector<unsigned> >
          sorted_index_groups_owner(index_groups.size());
        af::ref<std::vector<unsigned> >
          sorted_index_groups = sorted_index_groups_owner.ref();
        for(std::size_t i_cluster=0;i_cluster<perm.size();i_cluster++) {
          // swap is used to avoid deep copies of index_groups
          sorted_index_groups[i_cluster].swap(index_groups[perm[i_cluster]]);
        }
        index_groups = sorted_index_groups_owner;
        return *this;
      }

      //! In-place sorting of index_groups.
      asu_clusters&
      sort_indices_in_each_group()
      {
        af::ref<std::vector<unsigned> > index_groups_ref = index_groups.ref();
        for(std::size_t
              i_cluster=0;
              i_cluster<index_groups_ref.size();
              i_cluster++){
          std::sort(
            index_groups_ref[i_cluster].begin(),
            index_groups_ref[i_cluster].end());
        }
        return *this;
      }

      /*! \brief Indices with respect to the asu_mappings of the
          pair_asu_table as passed to the constructor.
       */
      af::shared<std::vector<unsigned> > index_groups;
  };

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_ASU_CLUSTERS_H
