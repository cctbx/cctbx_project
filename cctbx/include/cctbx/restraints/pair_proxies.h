#ifndef CCTBX_RESTRAINTS_PAIR_PROXIES_H
#define CCTBX_RESTRAINTS_PAIR_PROXIES_H

#include <cctbx/restraints/bond.h>
#include <cctbx/restraints/repulsion.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace restraints {

  class pair_proxies
  {
    public:
      pair_proxies() {}

      template <typename ScattererType>
      pair_proxies(
        af::const_ref<ScattererType> const& scatterers,
        restraints::bond_params_table const& bond_params_table,
        restraints::repulsion_distance_table const& repulsion_distance_table,
        std::vector<crystal::pair_asu_table<> > const& shell_asu_tables,
        af::const_ref<double> const& shell_distance_cutoffs,
        double nonbonded_distance_cutoff,
        double nonbonded_buffer,
        double vdw_1_4_factor)
      {
        CCTBX_ASSERT(shell_asu_tables.size() > 0);
        CCTBX_ASSERT(shell_distance_cutoffs.size() == shell_asu_tables.size());
        double distance_cutoff = std::max(
          af::max(shell_distance_cutoffs),
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
            bond_asu_proxies.push_back(
              make_bond_asu_proxy(bond_params_table, pair));
          }
          else if (shell_asu_tables.size() > 1
                   && shell_asu_tables[1].contains(pair)) {
            continue;
          }
          else if (shell_asu_tables.size() > 2
                   && shell_asu_tables[2].contains(pair)) {
            repulsion_asu_proxies.push_back(make_repulsion_asu_proxy(
              scatterers, repulsion_distance_table, pair, vdw_1_4_factor));
          }
          else if (pair.dist_sq <= nonbonded_distance_cutoff_sq) {
            repulsion_asu_proxies.push_back(make_repulsion_asu_proxy(
              scatterers, repulsion_distance_table, pair));
          }
        }
      }

      static
      bond_asu_proxy
      make_bond_asu_proxy(
        restraints::bond_params_table const& bond_params_table,
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

      template <typename ScattererType>
      static
      repulsion_asu_proxy
      make_repulsion_asu_proxy(
        af::const_ref<ScattererType> const& scatterers,
        restraints::repulsion_distance_table const& repulsion_distance_table,
        direct_space_asu::asu_mapping_index_pair const& pair,
        double vdw_factor=1)
      {
        std::string const& sc_type_i = scatterers[pair.i_seq].scattering_type;
        std::string const& sc_type_j = scatterers[pair.j_seq].scattering_type;
        repulsion_distance_table::const_iterator
          distance_dict = repulsion_distance_table.find(sc_type_i);
        static const char* error_msg =
         "Unknown scattering type pair (incomplete repulsion_distance_table).";
        if (distance_dict == repulsion_distance_table.end()) {
          throw error(error_msg);
        }
        repulsion_distance_dict::const_iterator
          dict_entry = distance_dict->second.find(sc_type_j);
        return repulsion_asu_proxy(pair, dict_entry->second*vdw_factor);
      }

      af::shared<restraints::bond_asu_proxy> bond_asu_proxies;
      af::shared<restraints::repulsion_asu_proxy> repulsion_asu_proxies;
  };

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_PAIR_PROXIES_H
