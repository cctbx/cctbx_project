#ifndef CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H
#define CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H

#include <cctbx/crystal/direct_space_asu.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace crystal {

//! Coordination sequence algorithms.
namespace coordination_sequences {

  /*! \brief Grouping of symmetry operations characterizing a node in a
      graph of bonded sites.
   */
  struct node
  {
    //! Default constructor. Some data members are not initialized!
    node() {}

    //! Initialization of rt_mx and rt_mx_unique.
    node(
      sgtbx::rt_mx const& special_op,
      sgtbx::rt_mx const& rt_mx_)
    :
      rt_mx(rt_mx_),
      rt_mx_unique(rt_mx_.multiply(special_op))
    {}

    //! Matrix as passed to the constructor.
    sgtbx::rt_mx rt_mx;
    //! rt_mx * special_op(i_seq).
    sgtbx::rt_mx rt_mx_unique;
  };

  //! Buffer for nodes in three shells.
  struct three_shells
  {
    //! Initialization of three empty shells.
    three_shells()
    {
      prev   = &all[0];
      middle = &all[1];
      next   = &all[2];
    }

    //! To be called when starting with a new pivot site.
    /*! Clearing of the three shells. The next shell is initialized
        with the identity matrix (for the pivot site).
     */
    void
    clear(
      sgtbx::rt_mx const& special_op,
      unsigned i_seq_pivot)
    {
      prev->clear();
      middle->clear();
      next->clear();
      (*next)[i_seq_pivot].push_back(node(special_op, sgtbx::rt_mx(1,1)));
    }

    //! To be called when starting with a new shell.
    void
    shift()
    {
      std::map<unsigned, std::vector<node> >* tmp = prev;
      prev = middle;
      middle = next;
      next = tmp;
      next->clear();
    }

    //! Number of nodes in next shell.
    unsigned
    count_next() const
    {
      unsigned result = 0;
      for(std::map<unsigned, std::vector<node> >::const_iterator
            item_n=next->begin();
            item_n!=next->end();
            item_n++) {
        result += item_n->second.size();
      }
      return result;
    }

    //! Search for node in all three shells.
    /*! Returns true if the test_node was found, false otherwise.
     */
    bool
    find_node(unsigned j_seq, node const& test_node) const
    {
      for(unsigned i_shell=0;i_shell<3;i_shell++) {
        std::map<unsigned, std::vector<node> > const& shell = all[i_shell];
        std::map<unsigned, std::vector<node> >::const_iterator
          item = shell.find(j_seq);
        if (item == shell.end()) continue;
        std::vector<node> const& list_nodes = item->second;
        for(std::vector<node>::const_iterator
              list_node=list_nodes.begin();
              list_node!=list_nodes.end();
              list_node++) {
          if (list_node->rt_mx_unique == test_node.rt_mx_unique) {
            return true;
          }
        }
      }
      return false;
    }

    //! Memory for all three shells.
    af::tiny<std::map<unsigned, std::vector<node> >, 3> all;
    //! Pointer to memory for prev shell.
    std::map<unsigned, std::vector<node> >* prev;
    //! Pointer to memory for middle shell.
    std::map<unsigned, std::vector<node> >* middle;
    //! Pointer to memory for next shell.
    std::map<unsigned, std::vector<node> >* next;
  };

  //! Generic coordination sequence algorithm.
  template <typename Actions>
  struct core_asu : Actions
  {
    //! Execution of the coordination sequence algorithm.
    core_asu(
      crystal::pair_asu_table<> const& pair_asu_table,
      unsigned max_shell)
    :
      Actions(pair_asu_table, max_shell)
    {
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *pair_asu_table.asu_mappings().get();
      af::const_ref<pair_asu_dict>
        asu_table_ref = pair_asu_table.table().const_ref();
      CCTBX_ASSERT(asu_mappings.mappings_const_ref().size()
                == asu_table_ref.size());
      unsigned n_seq = asu_table_ref.size();
      three_shells shells;
      for(this->i_seq_pivot=0; this->i_seq_pivot<n_seq; this->i_seq_pivot++) {
        shells.clear(
          asu_mappings.special_op(this->i_seq_pivot),
          this->i_seq_pivot);
        for(this->i_shell_minus_1=0;
            this->i_shell_minus_1<max_shell;
            this->i_shell_minus_1++) {
          shells.shift();
          for(std::map<unsigned, std::vector<node> >::const_iterator
                items_m=shells.middle->begin();
                items_m!=shells.middle->end();
                items_m++) {
            unsigned i_seq_m = items_m->first;
            std::vector<node> const& nodes_middle = items_m->second;
            sgtbx::rt_mx
              rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq_m, 0).inverse();
            for(unsigned i_node_m=0;i_node_m<nodes_middle.size();i_node_m++) {
              node node_m = nodes_middle[i_node_m];
              sgtbx::rt_mx rt_mx_ni = node_m.rt_mx.multiply(rt_mx_i_inv);
              pair_asu_dict::const_iterator
                pair_asu_dict_end = asu_table_ref[i_seq_m].end();
              for(pair_asu_dict::const_iterator
                    pair_asu_dict_i = asu_table_ref[i_seq_m].begin();
                    pair_asu_dict_i != pair_asu_dict_end;
                    pair_asu_dict_i++) {
                unsigned j_seq = pair_asu_dict_i->first;
                pair_asu_j_sym_groups const&
                  j_sym_groups=pair_asu_dict_i->second;
                for(unsigned i_group=0;
                             i_group<j_sym_groups.size();
                             i_group++) {
                  pair_asu_j_sym_group j_sym_group = j_sym_groups[i_group];
                  pair_asu_j_sym_group::const_iterator
                    j_sym_group_end = j_sym_group.end();
                  for(pair_asu_j_sym_group::const_iterator
                        j_sym_group_i = j_sym_group.begin();
                        j_sym_group_i != j_sym_group_end;
                        j_sym_group_i++) {
                    unsigned j_sym = *j_sym_group_i;
                    sgtbx::rt_mx
                      rt_mx_j = asu_mappings.get_rt_mx(j_seq, j_sym);
                    node new_node(
                      asu_mappings.special_op(j_seq),
                      rt_mx_ni.multiply(rt_mx_j));
                    if (!shells.find_node(j_seq, new_node)) {
                      (*shells.next)[j_seq].push_back(new_node);
                    }
                  }
                }
              }
            }
          }
          this->shell_complete(shells);
        }
        this->pivot_complete();
      }
    }
  };

  //! Generic coordination sequence algorithm.
  template <typename Actions>
  struct core_sym : Actions
  {
    //! Execution of the coordination sequence algorithm.
    core_sym(
      crystal::pair_sym_table const& full_pair_sym_table,
      sgtbx::site_symmetry_table const& site_symmetry_table,
      unsigned max_shell)
    :
      Actions(full_pair_sym_table, max_shell)
    {
      CCTBX_ASSERT(
        full_pair_sym_table.size() == site_symmetry_table.indices().size());
      unsigned n_seq = full_pair_sym_table.size();
      three_shells shells;
      for(this->i_seq_pivot=0; this->i_seq_pivot<n_seq; this->i_seq_pivot++) {
        shells.clear(
          site_symmetry_table.get(this->i_seq_pivot).special_op(),
          this->i_seq_pivot);
        for(this->i_shell_minus_1=0;
            this->i_shell_minus_1<max_shell;
            this->i_shell_minus_1++) {
          shells.shift();
          for(std::map<unsigned, std::vector<node> >::const_iterator
                items_m=shells.middle->begin();
                items_m!=shells.middle->end();
                items_m++) {
            unsigned i_seq_m = items_m->first;
            std::vector<node> const& nodes_middle = items_m->second;
            unsigned nodes_middle_size = nodes_middle.size();
            for(unsigned i_node_m=0;i_node_m<nodes_middle_size;i_node_m++) {
              node node_m = nodes_middle[i_node_m];
              pair_sym_dict::const_iterator
                pair_sym_dict_end = full_pair_sym_table[i_seq_m].end();
              for(pair_sym_dict::const_iterator
                    pair_sym_dict_i = full_pair_sym_table[i_seq_m].begin();
                    pair_sym_dict_i != pair_sym_dict_end;
                    pair_sym_dict_i++) {
                unsigned j_seq = pair_sym_dict_i->first;
                pair_sym_ops const& j_sym_ops = pair_sym_dict_i->second;
                for(unsigned i_j_sym_op=0;
                             i_j_sym_op<j_sym_ops.size();
                             i_j_sym_op++) {
                  node new_node(
                    site_symmetry_table.get(j_seq).special_op(),
                    node_m.rt_mx.multiply(j_sym_ops[i_j_sym_op]));
                  if (!shells.find_node(j_seq, new_node)) {
                    (*shells.next)[j_seq].push_back(new_node);
                  }
                }
              }
            }
          }
          this->shell_complete(shells);
        }
        this->pivot_complete();
      }
    }
  };

  //! Actions for simple counting.
  struct term_table_actions_asu
  {
    //! Called at start of core_asu<term_table_actions_asu>.
    term_table_actions_asu(
      crystal::pair_asu_table<> const& pair_asu_table,
      unsigned /*max_shell*/)
    :
      terms(0)
    {
      term_table.reserve(pair_asu_table.table().size());
    }

    //! Called when the next shell is complete.
    void
    shell_complete(three_shells const& shells)
    {
      if (terms == 0) {
        term_table.push_back(std::vector<unsigned>());
        terms = &term_table.back();
        terms->push_back(1);
      }
      terms->push_back(shells.count_next());
    }

    //! Called when a pivot site is completed.
    void
    pivot_complete() { terms = 0; }

    //! Index of current pivot site.
    unsigned i_seq_pivot;
    //! Index of current shell - 1.
    unsigned i_shell_minus_1;
    //! Final table of terms.
    af::shared<std::vector<unsigned> > term_table;
    //! Pointer to term_table[i_seq_pivot].
    std::vector<unsigned>* terms;
  };

  //! Actions for simple counting.
  struct term_table_actions_sym
  {
    //! Called at start of core_sym<term_table_actions_sym>.
    term_table_actions_sym(
      crystal::pair_sym_table const& full_pair_sym_table,
      unsigned /*max_shell*/)
    :
      terms(0)
    {
      term_table.reserve(full_pair_sym_table.size());
    }

    //! Called when the next shell is complete.
    void
    shell_complete(three_shells const& shells)
    {
      if (terms == 0) {
        term_table.push_back(std::vector<unsigned>());
        terms = &term_table.back();
        terms->push_back(1);
      }
      terms->push_back(shells.count_next());
    }

    //! Called when a pivot site is completed.
    void
    pivot_complete() { terms = 0; }

    //! Index of current pivot site.
    unsigned i_seq_pivot;
    //! Index of current shell - 1.
    unsigned i_shell_minus_1;
    //! Final table of terms.
    af::shared<std::vector<unsigned> > term_table;
    //! Pointer to term_table[i_seq_pivot].
    std::vector<unsigned>* terms;
  };

  //! Friendly interface to core_asu<term_table_actions_asu>.
  af::shared<std::vector<unsigned> >
  simple(
    crystal::pair_asu_table<> const& pair_asu_table,
    unsigned max_shell)
  {
    return core_asu<term_table_actions_asu>(
      pair_asu_table, max_shell).term_table;
  }

  //! Friendly interface to core_asu<term_table_actions_asu>.
  af::shared<std::vector<unsigned> >
  simple_sym(
    crystal::pair_sym_table const& full_pair_sym_table,
    sgtbx::site_symmetry_table const& site_symmetry_table,
    unsigned max_shell)
  {
    return core_sym<term_table_actions_sym>(
      full_pair_sym_table, site_symmetry_table, max_shell).term_table;
  }

  //! Actions for the generation of higher-level (nonbonded) interactions.
  struct shell_asu_tables_actions
  {
    //! Called at start of core_asu<shell_asu_tables_actions>.
    shell_asu_tables_actions(
      crystal::pair_asu_table<> const& pair_asu_table,
      unsigned max_shell)
    {
      shell_asu_tables.reserve(max_shell);
      if (max_shell > 0) {
        shell_asu_tables.push_back(pair_asu_table);
        for(i_shell_minus_1=1;
            i_shell_minus_1<max_shell;
            i_shell_minus_1++) {
          shell_asu_tables.push_back(
            crystal::pair_asu_table<>(pair_asu_table.asu_mappings()));
        }
      }
    }

    //! Called when the next shell is complete.
    void
    shell_complete(three_shells const& shells)
    {
      if (i_shell_minus_1 == 0) return;
      pair_asu_table<>& shell_asu_table = shell_asu_tables[i_shell_minus_1];
      direct_space_asu::asu_mappings<>*
        asu_mappings = shell_asu_table.asu_mappings().get();
      sgtbx::rt_mx rt_mx_i = asu_mappings->get_rt_mx(i_seq_pivot, 0);
      for(std::map<unsigned, std::vector<node> >::const_iterator
            items_n=shells.next->begin();
            items_n!=shells.next->end();
            items_n++) {
        unsigned i_seq_node = items_n->first;
        std::vector<node> const& nodes = items_n->second;
        for(unsigned i_node=0;i_node<nodes.size();i_node++) {
          node const& node_ = nodes[i_node];
          sgtbx::rt_mx rt_mx_j = rt_mx_i.multiply(node_.rt_mx);
          int j_sym = asu_mappings->find_i_sym(i_seq_node, rt_mx_j);
          if (j_sym < 0) continue;
          shell_asu_table.process_pair(
            i_seq_pivot,
            i_seq_node,
            /* rt_mx_ji = */ node_.rt_mx,
            rt_mx_i,
            j_sym);
        }
      }
    }

    //! Called when a pivot site is completed.
    void
    pivot_complete() {}

    //! Index of current pivot site.
    unsigned i_seq_pivot;
    //! Index of current shell - 1.
    unsigned i_shell_minus_1;
    //! Final array of pair_asu_table instances.
    std::vector<pair_asu_table<> > shell_asu_tables;
  };

  //! Friendly interface to core_asu<shell_asu_tables_actions>.
  std::vector<pair_asu_table<> >
  shell_asu_tables(
    crystal::pair_asu_table<> const& pair_asu_table,
    unsigned max_shell)
  {
    return core_asu<shell_asu_tables_actions>(pair_asu_table, max_shell)
      .shell_asu_tables;
  }

}}} // namespace cctbx::crystal::coordination_sequences

#endif // CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H
