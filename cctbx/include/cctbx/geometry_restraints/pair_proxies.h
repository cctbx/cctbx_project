#ifndef CCTBX_GEOMETRY_RESTRAINTS_PAIR_PROXIES_H
#define CCTBX_GEOMETRY_RESTRAINTS_PAIR_PROXIES_H

#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry_restraints/nonbonded.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace geometry_restraints {

  //! Adds all pairs defined by bond_simple_proxies[i].i_seqs.
  void
  add_pairs(
    crystal::pair_asu_table<>& pair_asu_table,
    af::const_ref<bond_simple_proxy> const& bond_simple_proxies)
  {
    for(unsigned i=0;i<bond_simple_proxies.size();i++) {
      pair_asu_table.add_pair(bond_simple_proxies[i].i_seqs);
    }
  }

  //! Generation of sorted bond and respulsion proxies.
  class pair_proxies
  {
    public:
      //! Default constructor. Some data members are not initialized!
      pair_proxies() {}

      //! Generation of bond proxies only.
      pair_proxies(
        af::const_ref<bond_params_dict> const& bond_params_table)
      :
        n_bonded(0),
        n_1_3(0),
        n_1_4(0),
        n_nonbonded(0),
        n_unknown_nonbonded_type_pairs(0)
      {
        for(unsigned i_seq=0;i_seq<bond_params_table.size();i_seq++) {
          for(bond_params_dict::const_iterator
                dict_i=bond_params_table[i_seq].begin();
                dict_i!=bond_params_table[i_seq].end();
                dict_i++) {
            bond_proxies.process(bond_simple_proxy(
              af::tiny<unsigned, 2>(i_seq, dict_i->first),
              dict_i->second));
          }
        }
      }

      //! Generation of bond and nonbonded proxies.
      /*! Bonds (1-2) interactions are defined by shell_asu_tables[0]
          and are used at the same time as nonbonded exclusions.

          1-3 interactions are defined by shell_asu_tables[1]
          and are used as nonbonded exclusions.

          1-4 interactions are defined by shell_asu_tables[2]
          and are used to attenuate nonbonded energies according to
          nonbonded_params.
       */
      pair_proxies(
        af::const_ref<std::size_t> const& model_indices,
        af::const_ref<std::size_t> const& conformer_indices,
        geometry_restraints::nonbonded_params const& nonbonded_params,
        af::const_ref<std::string> const& nonbonded_types,
        af::const_ref<bond_params_dict> const& bond_params_table,
        std::vector<crystal::pair_asu_table<> > const& shell_asu_tables,
        double bonded_distance_cutoff,
        double nonbonded_distance_cutoff,
        double nonbonded_buffer)
      :
        n_bonded(0),
        n_1_3(0),
        n_1_4(0),
        n_nonbonded(0),
        n_unknown_nonbonded_type_pairs(0)
      {
        CCTBX_ASSERT(model_indices.size() == 0
                  || model_indices.size() == nonbonded_types.size());
        CCTBX_ASSERT(conformer_indices.size() == 0
                  || conformer_indices.size() == nonbonded_types.size());
        CCTBX_ASSERT(bond_params_table.size() == nonbonded_types.size());
        CCTBX_ASSERT(shell_asu_tables.size() > 0);
        for(unsigned i=0; i<shell_asu_tables.size(); i++) {
          CCTBX_ASSERT(shell_asu_tables[i].table().size()
                    == bond_params_table.size());
        }
        for(unsigned i=1; i<shell_asu_tables.size(); i++) {
          CCTBX_ASSERT(shell_asu_tables[i].asu_mappings().get()
                    == shell_asu_tables[0].asu_mappings().get());
        }
        bond_proxies = bond_sorted_asu_proxies(
          shell_asu_tables[0].asu_mappings());
        nonbonded_proxies = nonbonded_sorted_asu_proxies(
          shell_asu_tables[0].asu_mappings());
        double distance_cutoff = std::max(
          bonded_distance_cutoff,
          nonbonded_distance_cutoff+nonbonded_buffer);
        bool minimal = false;
        crystal::neighbors::fast_pair_generator<> pair_generator(
          shell_asu_tables[0].asu_mappings(),
          distance_cutoff,
          minimal);
        double nonbonded_distance_cutoff_sq = nonbonded_distance_cutoff
                                            * nonbonded_distance_cutoff;
        while (!pair_generator.at_end()) {
          direct_space_asu::asu_mapping_index_pair_and_diff<>
            pair = pair_generator.next();
          if (shell_asu_tables[0].contains(pair)) {
            bond_proxies.process(
              make_bond_asu_proxy(bond_params_table, pair));
            n_bonded++;
          }
          else if (shell_asu_tables.size() > 1
                   && shell_asu_tables[1].contains(pair)) {
            n_1_3++;
          }
          else if ((model_indices.size() == 0
                    || model_indices[pair.i_seq] == model_indices[pair.j_seq])
                   && (conformer_indices.size() == 0
                       || conformer_indices[pair.i_seq] == 0
                       || conformer_indices[pair.j_seq] == 0
                       || conformer_indices[pair.i_seq]
                       == conformer_indices[pair.j_seq])) {
            if (shell_asu_tables.size() > 2
                && shell_asu_tables[2].contains(pair)) {
              nonbonded_proxies.process(make_nonbonded_asu_proxy(
                nonbonded_params, nonbonded_types, pair, true));
              n_1_4++;
            }
            else if (pair.dist_sq <= nonbonded_distance_cutoff_sq) {
              nonbonded_proxies.process(make_nonbonded_asu_proxy(
                nonbonded_params, nonbonded_types, pair, false));
              n_nonbonded++;
            }
          }
        }
      }

      //! For internal use only.
      static
      bond_asu_proxy
      make_bond_asu_proxy(
        af::const_ref<bond_params_dict> const& bond_params_table,
        direct_space_asu::asu_mapping_index_pair const& pair)
      {
        unsigned i, j;
        if (pair.i_seq <= pair.j_seq) {
          i = pair.i_seq;
          j = pair.j_seq;
        }
        else {
          i = pair.j_seq;
          j = pair.i_seq;
        }
        bond_params_dict const& params_dict = bond_params_table[i];
        bond_params_dict::const_iterator params = params_dict.find(j);
        if (params == params_dict.end()) {
          throw error(
            "Unknown bond parameters (incomplete bond_params_table).");
        }
        return bond_asu_proxy(pair, params->second);
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

      //! Final bond proxies.
      bond_sorted_asu_proxies bond_proxies;
      //! Final nonbonded proxies.
      nonbonded_sorted_asu_proxies nonbonded_proxies;
      //! Total number of bond proxies (simple + asu).
      unsigned n_bonded;
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

#endif // CCTBX_GEOMETRY_RESTRAINTS_PAIR_PROXIES_H
