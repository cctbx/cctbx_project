#ifndef CCTBX_GEOMETRY_RESTRAINTS_PAIR_PROXIES_H
#define CCTBX_GEOMETRY_RESTRAINTS_PAIR_PROXIES_H

#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry_restraints/repulsion.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace geometry_restraints {

  void
  add_pairs(
    crystal::pair_asu_table<>& pair_asu_table,
    af::const_ref<bond_simple_proxy> const& bond_simple_proxies)
  {
    for(unsigned i=0;i<bond_simple_proxies.size();i++) {
      pair_asu_table.add_pair(bond_simple_proxies[i].i_seqs);
    }
  }

  class pair_proxies
  {
    public:
      pair_proxies() {}

      pair_proxies(
        af::const_ref<bond_params_dict> const& bond_params_table)
      :
        n_bonded(0),
        n_1_3(0),
        n_1_4(0),
        n_nonbonded(0),
        n_unknown_repulsion_type_pairs(0)
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

      pair_proxies(
        geometry_restraints::repulsion_params const& repulsion_params,
        af::const_ref<std::string> const& repulsion_types,
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
        n_unknown_repulsion_type_pairs(0)
      {
        CCTBX_ASSERT(repulsion_types.size() == bond_params_table.size());
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
        repulsion_proxies = repulsion_sorted_asu_proxies(
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
          else if (shell_asu_tables.size() > 2
                   && shell_asu_tables[2].contains(pair)) {
            repulsion_proxies.process(make_repulsion_asu_proxy(
              repulsion_params, repulsion_types, pair, true));
            n_1_4++;
          }
          else if (pair.dist_sq <= nonbonded_distance_cutoff_sq) {
            repulsion_proxies.process(make_repulsion_asu_proxy(
              repulsion_params, repulsion_types, pair, false));
            n_nonbonded++;
          }
        }
      }

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

      repulsion_asu_proxy
      make_repulsion_asu_proxy(
        geometry_restraints::repulsion_params const& repulsion_params,
        af::const_ref<std::string> const& repulsion_types,
        direct_space_asu::asu_mapping_index_pair const& pair,
        bool is_1_4_interaction)
      {
        std::string const& rep_type_i = repulsion_types[pair.i_seq];
        std::string const& rep_type_j = repulsion_types[pair.j_seq];
        repulsion_distance_table::const_iterator
          distance_dict = repulsion_params.distance_table.find(rep_type_i);
        if (distance_dict != repulsion_params.distance_table.end()) {
          repulsion_distance_dict::const_iterator
            dict_entry = distance_dict->second.find(rep_type_j);
          if (dict_entry != distance_dict->second.end()) {
            return repulsion_asu_proxy(pair, get_repulsion_distance(
              repulsion_params, dict_entry->second, is_1_4_interaction));
          }
        }
        distance_dict = repulsion_params.distance_table.find(rep_type_j);
        if (distance_dict != repulsion_params.distance_table.end()) {
          repulsion_distance_dict::const_iterator
            dict_entry = distance_dict->second.find(rep_type_i);
          if (dict_entry != distance_dict->second.end()) {
            return repulsion_asu_proxy(pair, get_repulsion_distance(
              repulsion_params, dict_entry->second, is_1_4_interaction));
          }
        }
        geometry_restraints::repulsion_radius_table::const_iterator
          radius_i = repulsion_params.radius_table.find(rep_type_i);
        if (radius_i != repulsion_params.radius_table.end()) {
          geometry_restraints::repulsion_radius_table::const_iterator
            radius_j = repulsion_params.radius_table.find(rep_type_j);
          if (radius_j != repulsion_params.radius_table.end()) {
            return repulsion_asu_proxy(pair, get_repulsion_distance(
              repulsion_params,
              radius_i->second + radius_j->second,
              is_1_4_interaction));
          }
        }
        if (repulsion_params.default_distance > 0) {
          n_unknown_repulsion_type_pairs++;
          return repulsion_asu_proxy(pair, get_repulsion_distance(
            repulsion_params,
            repulsion_params.default_distance,
            is_1_4_interaction));
        }
        throw error(
         "Unknown repulsion type pair (incomplete repulsion_distance_table): "
         + rep_type_i + " - " + rep_type_j);
      }

      double
      get_repulsion_distance(
        geometry_restraints::repulsion_params const& repulsion_params,
        double distance,
        bool is_1_4_interaction)
      {
        if (is_1_4_interaction) {
          distance *= repulsion_params.factor_1_4_interactions;
          distance -= repulsion_params.const_shrink_1_4_interactions;
        }
        return std::max(repulsion_params.minimum_distance, distance);
      }

      bond_sorted_asu_proxies bond_proxies;
      repulsion_sorted_asu_proxies repulsion_proxies;
      unsigned n_bonded;
      unsigned n_1_3;
      unsigned n_1_4;
      unsigned n_nonbonded;
      unsigned n_unknown_repulsion_type_pairs;
  };

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_PAIR_PROXIES_H
