#ifndef CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H
#define CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H

#include <cctbx/crystal/direct_space_asu.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace crystal {

//! Coordination sequence algorithms.
namespace coordination_sequences {

  struct node
  {
    node() {}

    node(
      direct_space_asu::asu_mappings<> const& asu_mappings,
      unsigned i_seq,
      sgtbx::rt_mx const& rt_mx_)
    :
      rt_mx(rt_mx_),
      rt_mx_unique(rt_mx_.multiply(asu_mappings.special_op(i_seq)))
    {}

    sgtbx::rt_mx rt_mx;
    sgtbx::rt_mx rt_mx_unique;
  };

  struct three_shells
  {
    three_shells()
    {
      prev   = &all[0];
      middle = &all[1];
      next   = &all[2];
    }

    void
    clear(
      direct_space_asu::asu_mappings<> const& asu_mappings,
      unsigned i_seq_pivot)
    {
      prev->clear();
      middle->clear();
      next->clear();
      (*next)[i_seq_pivot].push_back(
        node(asu_mappings, i_seq_pivot, sgtbx::rt_mx(1,1)));
    }

    void
    shift()
    {
      std::map<unsigned, std::vector<node> >* tmp = prev;
      prev = middle;
      middle = next;
      next = tmp;
      next->clear();
    }

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

    bool
    find_node(unsigned j_seq, node const& test_node) const
    {
      for(unsigned i_shell=0;i_shell<3;i_shell++) {
        std::map<unsigned, std::vector<node> > const& shell = all[i_shell];
        std::map<unsigned, std::vector<node> >::const_iterator
          item = shell.find(j_seq);
        if (item == shell.end()) continue;
        for(const node*
              list_node=&*item->second.begin();
              list_node!=&*item->second.end();
              list_node++) {
          if (list_node->rt_mx_unique == test_node.rt_mx_unique) {
            return true;
          }
        }
      }
      return false;
    }

    af::tiny<std::map<unsigned, std::vector<node> >, 3> all;
    std::map<unsigned, std::vector<node> >* prev;
    std::map<unsigned, std::vector<node> >* middle;
    std::map<unsigned, std::vector<node> >* next;
  };

  template <typename Actions>
  struct core : Actions
  {
    core(
      crystal::pair_asu_table<> const& pair_asu_table,
      unsigned n_shells)
    :
      Actions(pair_asu_table, n_shells)
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
        sgtbx::rt_mx rt_mx_pivot = asu_mappings.get_rt_mx(this->i_seq_pivot,0);
        shells.clear(asu_mappings, this->i_seq_pivot);
        for(this->i_shell_minus_1=0;
            this->i_shell_minus_1<n_shells;
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
                      asu_mappings, j_seq, rt_mx_ni.multiply(rt_mx_j));
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

  struct term_table_actions
  {
    term_table_actions(
      crystal::pair_asu_table<> const& pair_asu_table,
      unsigned n_shells)
    :
      terms(0)
    {
      term_table.reserve(pair_asu_table.table().size());
    }

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

    void
    pivot_complete() { terms = 0; }

    unsigned i_seq_pivot;
    unsigned i_shell_minus_1;
    af::shared<std::vector<unsigned> > term_table;
    std::vector<unsigned>* terms;
  };

  af::shared<std::vector<unsigned> >
  simple(
    crystal::pair_asu_table<> const& pair_asu_table,
    unsigned n_shells)
  {
    return core<term_table_actions>(pair_asu_table, n_shells).term_table;
  }

  struct shell_asu_tables_actions
  {
    shell_asu_tables_actions(
      crystal::pair_asu_table<> const& pair_asu_table,
      unsigned n_shells)
    {
      shell_asu_tables.reserve(n_shells);
      if (n_shells > 0) {
        shell_asu_tables.push_back(pair_asu_table);
        for(i_shell_minus_1=1;
            i_shell_minus_1<n_shells;
            i_shell_minus_1++) {
          shell_asu_tables.push_back(
            crystal::pair_asu_table<>(pair_asu_table.asu_mappings()));
        }
      }
    }

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

    void
    pivot_complete() {}

    unsigned i_seq_pivot;
    unsigned i_shell_minus_1;
    std::vector<pair_asu_table<> > shell_asu_tables;
  };

  std::vector<pair_asu_table<> >
  shell_asu_tables(
    crystal::pair_asu_table<> const& pair_asu_table,
    unsigned n_shells)
  {
    return core<shell_asu_tables_actions>(pair_asu_table, n_shells)
      .shell_asu_tables;
  }

}}} // namespace cctbx::crystal::coordination_sequences

#endif // CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H
