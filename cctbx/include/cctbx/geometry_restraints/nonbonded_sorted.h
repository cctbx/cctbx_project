#ifndef CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_SORTED_H
#define CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_SORTED_H

#include <cctbx/geometry_restraints/nonbonded.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace geometry_restraints {

  //! Generation of sorted nonbonded proxies.
  class nonbonded_sorted_asu_proxies : public nonbonded_sorted_asu_proxies_base
  {
    public:
      //! Default constructor. Some data members are not initialized!
      nonbonded_sorted_asu_proxies() {}

      //! Initialization with asu_mappings.
      nonbonded_sorted_asu_proxies(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<> > const& asu_mappings)
      :
        nonbonded_sorted_asu_proxies_base(asu_mappings),
        n_1_3(0),
        n_1_4(0),
        n_nonbonded(0),
        n_unknown_nonbonded_type_pairs(0)
      {}

      //! Generation of nonbonded proxies.
      /*! Bonds (1-2) interactions are defined by shell_asu_tables[0]
          and are used at the same time as nonbonded exclusions.

          1-3 interactions are defined by shell_asu_tables[1]
          and are used as nonbonded exclusions.

          1-4 interactions are defined by shell_asu_tables[2]
          and are used to attenuate nonbonded energies according to
          nonbonded_params.
       */
      nonbonded_sorted_asu_proxies(
        af::const_ref<std::size_t> const& model_indices,
        af::const_ref<std::size_t> const& conformer_indices,
        geometry_restraints::nonbonded_params const& nonbonded_params,
        af::const_ref<std::string> const& nonbonded_types,
        double nonbonded_distance_cutoff_plus_buffer,
        std::vector<crystal::pair_asu_table<> > const& shell_asu_tables)
      :
        nonbonded_sorted_asu_proxies_base(shell_asu_tables[0].asu_mappings()),
        n_1_3(0),
        n_1_4(0),
        n_nonbonded(0),
        n_unknown_nonbonded_type_pairs(0)
      {
        CCTBX_ASSERT(model_indices.size() == 0
                  || model_indices.size() == nonbonded_types.size());
        CCTBX_ASSERT(conformer_indices.size() == 0
                  || conformer_indices.size() == nonbonded_types.size());
        CCTBX_ASSERT(shell_asu_tables.size() > 0);
        for(unsigned i=0; i<shell_asu_tables.size(); i++) {
          CCTBX_ASSERT(shell_asu_tables[i].table().size()
                    == nonbonded_types.size());
        }
        unsigned shell_asu_tables_size = shell_asu_tables.size();
        for(unsigned i=1; i<shell_asu_tables_size; i++) {
          CCTBX_ASSERT(shell_asu_tables[i].asu_mappings().get()
                    == shell_asu_tables[0].asu_mappings().get());
        }
        bool minimal = false;
        crystal::neighbors::fast_pair_generator<> pair_generator(
          shell_asu_tables[0].asu_mappings(),
          nonbonded_distance_cutoff_plus_buffer,
          minimal);
        while (!pair_generator.at_end()) {
          direct_space_asu::asu_mapping_index_pair_and_diff<>
            pair = pair_generator.next();
          if (shell_asu_tables[0].contains(pair)) {
            continue;
          }
          if (   shell_asu_tables_size > 1
              && shell_asu_tables[1].contains(pair)) {
            n_1_3++;
            continue;
          }
          if (   model_indices.size() != 0
              && model_indices[pair.i_seq] != model_indices[pair.j_seq]) {
            continue;
          }
          if (   conformer_indices.size() != 0
              && conformer_indices[pair.i_seq] != 0
              && conformer_indices[pair.j_seq] != 0
              && conformer_indices[pair.i_seq]
              != conformer_indices[pair.j_seq]) {
            continue;
          }
          if (   shell_asu_tables_size > 2
              && shell_asu_tables[2].contains(pair)) {
            process(make_nonbonded_asu_proxy(
              nonbonded_params, nonbonded_types, pair, true));
            n_1_4++;
            continue;
          }
          process(make_nonbonded_asu_proxy(
            nonbonded_params, nonbonded_types, pair, false));
          n_nonbonded++;
        }
      }

      //! For internal use only.
      nonbonded_asu_proxy
      make_nonbonded_asu_proxy(
        geometry_restraints::nonbonded_params const& nonbonded_params,
        af::const_ref<std::string> const& nonbonded_types,
        direct_space_asu::asu_mapping_index_pair const& pair,
        bool is_1_4_interaction)
      {
        std::string const& rep_type_i = nonbonded_types[pair.i_seq];
        std::string const& rep_type_j = nonbonded_types[pair.j_seq];
        nonbonded_distance_table::const_iterator
          distance_dict = nonbonded_params.distance_table.find(rep_type_i);
        if (distance_dict != nonbonded_params.distance_table.end()) {
          nonbonded_distance_dict::const_iterator
            dict_entry = distance_dict->second.find(rep_type_j);
          if (dict_entry != distance_dict->second.end()) {
            return nonbonded_asu_proxy(pair, get_nonbonded_distance(
              nonbonded_params, dict_entry->second, is_1_4_interaction));
          }
        }
        distance_dict = nonbonded_params.distance_table.find(rep_type_j);
        if (distance_dict != nonbonded_params.distance_table.end()) {
          nonbonded_distance_dict::const_iterator
            dict_entry = distance_dict->second.find(rep_type_i);
          if (dict_entry != distance_dict->second.end()) {
            return nonbonded_asu_proxy(pair, get_nonbonded_distance(
              nonbonded_params, dict_entry->second, is_1_4_interaction));
          }
        }
        geometry_restraints::nonbonded_radius_table::const_iterator
          radius_i = nonbonded_params.radius_table.find(rep_type_i);
        if (radius_i != nonbonded_params.radius_table.end()) {
          geometry_restraints::nonbonded_radius_table::const_iterator
            radius_j = nonbonded_params.radius_table.find(rep_type_j);
          if (radius_j != nonbonded_params.radius_table.end()) {
            return nonbonded_asu_proxy(pair, get_nonbonded_distance(
              nonbonded_params,
              radius_i->second + radius_j->second,
              is_1_4_interaction));
          }
        }
        if (nonbonded_params.default_distance > 0) {
          n_unknown_nonbonded_type_pairs++;
          return nonbonded_asu_proxy(pair, get_nonbonded_distance(
            nonbonded_params,
            nonbonded_params.default_distance,
            is_1_4_interaction));
        }
        throw error(
         "Unknown nonbonded type pair (incomplete nonbonded_distance_table): "
         + rep_type_i + " - " + rep_type_j);
      }

      //! For internal use only.
      double
      get_nonbonded_distance(
        geometry_restraints::nonbonded_params const& nonbonded_params,
        double distance,
        bool is_1_4_interaction)
      {
        if (is_1_4_interaction) {
          distance *= nonbonded_params.factor_1_4_interactions;
          distance -= nonbonded_params.const_shrink_1_4_interactions;
        }
        return std::max(nonbonded_params.minimum_distance, distance);
      }

      //! Number of 1-3 interactions excluded from nonbonded proxies.
      unsigned n_1_3;
      //! Number of attenuated 1-4 interactions.
      unsigned n_1_4;
      //! Total number of full strength (not attenuated) nonbonded terms.
      unsigned n_nonbonded;
      /*! \brief Number of unknown nonbonded type pairs (nonbonded proxies
          with default distance).
       */
      unsigned n_unknown_nonbonded_type_pairs;
  };

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_SORTED_H
