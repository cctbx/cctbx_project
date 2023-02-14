#ifndef CCTBX_GEOMETRY_RESTRAINTS_BOND_SORTED_H
#define CCTBX_GEOMETRY_RESTRAINTS_BOND_SORTED_H

#include <cctbx/geometry_restraints/bond.h>
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

  //! Generation of sorted bond proxies.
  class bond_sorted_asu_proxies : public bond_sorted_asu_proxies_base
  {
    public:
      //! Default constructor. Some data members are not initialized!
      bond_sorted_asu_proxies() {}

      //! Initialization with asu_mappings.
      bond_sorted_asu_proxies(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<> > const& asu_mappings)
      :
        bond_sorted_asu_proxies_base(asu_mappings)
      {}

      //! Initialization of simple proxies only.
      bond_sorted_asu_proxies(
        af::const_ref<bond_params_dict> const& bond_params_table)
      {
        for(unsigned i_seq=0;i_seq<bond_params_table.size();i_seq++) {
          for(bond_params_dict::const_iterator
                dict_i=bond_params_table[i_seq].begin();
                dict_i!=bond_params_table[i_seq].end();
                dict_i++) {
            process(bond_simple_proxy(
              af::tiny<unsigned, 2>(i_seq, dict_i->first),
              dict_i->second));
          }
        }
      }

      //! Initialization with bond_params_table and bond_asu_table.
      bond_sorted_asu_proxies(
        af::const_ref<bond_params_dict> const& bond_params_table,
        crystal::pair_asu_table<> const& bond_asu_table)
      :
        bond_sorted_asu_proxies_base(bond_asu_table.asu_mappings())
      {
        CCTBX_ASSERT(bond_asu_table.table().size()
                  == bond_params_table.size());
        af::const_ref<crystal::pair_asu_dict>
          asu_tab = bond_asu_table.table().const_ref();
        direct_space_asu::asu_mapping_index_pair pair;
        for(pair.i_seq=0;pair.i_seq<asu_tab.size();pair.i_seq++) {
          crystal::pair_asu_dict const& asu_dict = asu_tab[pair.i_seq];
          for(crystal::pair_asu_dict::const_iterator
                asu_dict_i=asu_dict.begin();
                asu_dict_i!=asu_dict.end();
                asu_dict_i++) {
            pair.j_seq = asu_dict_i->first;
            const bond_params_dict* params_dict;
            bond_params_dict::const_iterator params;
            if (pair.i_seq <= pair.j_seq ) {
              params_dict = &bond_params_table[pair.i_seq];
              params = params_dict->find(pair.j_seq);
            }
            else {
              params_dict = &bond_params_table[pair.j_seq];
              params = params_dict->find(pair.i_seq);
            }
            if (params == params_dict->end()) {
              unsigned i = pair.i_seq;
              unsigned j = pair.j_seq;
              if (i > j) std::swap(i, j);
              params_dict = &bond_params_table[j];
              params = params_dict->find(i);
              char buf[256];
              if (params == params_dict->end()) {
                std::snprintf(buf, sizeof(buf),
                  "Unknown bond parameters (incomplete bond_params_table):"
                  " i_seq=%d, j_seq=%d", i, j);
              }
              else {
                std::snprintf(buf, sizeof(buf),
                  "Improper bond_params_table (requirement i_seq <= j_seq):"
                  " i_seq=%d, j_seq=%d", j, i);
              }
              throw error(buf);
            }
            crystal::pair_asu_j_sym_groups const& jgs = asu_dict_i->second;
            for(unsigned ig=0;ig<jgs.size();ig++) {
              crystal::pair_asu_j_sym_group const& jg = jgs[ig];
              for(crystal::pair_asu_j_sym_group::const_iterator
                    jg_i=jg.begin();
                    jg_i!=jg.end();
                    jg_i++) {
                pair.j_sym = *jg_i;
                if (!pair.is_active()) continue;
                process(bond_asu_proxy(pair, params->second));
              }
            }
          }
        }
      }

      /*! \brief Initialization with pair_asu_table,
          distance_ideal=distance_model and weight=1.
       */
      bond_sorted_asu_proxies(
        crystal::pair_asu_table<> const& pair_asu_table)
      :
        bond_sorted_asu_proxies_base(pair_asu_table.asu_mappings())
      {
        af::const_ref<crystal::pair_asu_dict>
          asu_tab = pair_asu_table.table().const_ref();
        direct_space_asu::asu_mapping_index_pair pair;
        for(pair.i_seq=0;pair.i_seq<asu_tab.size();pair.i_seq++) {
          crystal::pair_asu_dict const& asu_dict = asu_tab[pair.i_seq];
          for(crystal::pair_asu_dict::const_iterator
                asu_dict_i=asu_dict.begin();
                asu_dict_i!=asu_dict.end();
                asu_dict_i++) {
            pair.j_seq = asu_dict_i->first;
            crystal::pair_asu_j_sym_groups const& jgs = asu_dict_i->second;
            for(unsigned ig=0;ig<jgs.size();ig++) {
              crystal::pair_asu_j_sym_group const& jg = jgs[ig];
              double distance_ideal = 0;
              for(crystal::pair_asu_j_sym_group::const_iterator
                    jg_i=jg.begin();
                    jg_i!=jg.end();
                    jg_i++) {
                pair.j_sym = *jg_i;
                distance_ideal += this->asu_mappings_->diff_vec(pair).length();
              }
              if (jg.size() != 0) distance_ideal /= jg.size();
              for(crystal::pair_asu_j_sym_group::const_iterator
                    jg_i=jg.begin();
                    jg_i!=jg.end();
                    jg_i++) {
                pair.j_sym = *jg_i;
                if (!pair.is_active()) continue;
                process(
                  bond_asu_proxy(pair, distance_ideal, /* weight */ 1.0));
              }
            }
          }
        }
      }
  };

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_BOND_SORTED_H
