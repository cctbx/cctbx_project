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
      unsigned i_seq_,
      sgtbx::rt_mx const& rt_mx_)
    :
      i_seq(i_seq_),
      rt_mx(rt_mx_)
    {
      rt_mx_unique = rt_mx_.multiply(asu_mappings.special_op(i_seq));
    }

    unsigned i_seq;
    sgtbx::rt_mx rt_mx;
    sgtbx::rt_mx rt_mx_unique;
  };

  bool
  find_node(
    node const& test_node,
    std::vector<node> const& node_list)
  {
    for(std::vector<node>::const_iterator
          list_node=node_list.begin();
          list_node!=node_list.end();
          list_node++) {
      if (list_node->rt_mx_unique == test_node.rt_mx_unique) {
        return true;
      }
    }
    return false;
  }

  struct three_shells
  {
    typedef std::vector<std::vector<node> > container;

    three_shells() {}

    three_shells(
      direct_space_asu::asu_mappings<> const& asu_mappings,
      unsigned i_seq_pivot)
    :
      a(asu_mappings.mappings().size()),
      b(asu_mappings.mappings().size()),
      c(asu_mappings.mappings().size())
    {
      prev = &a;
      middle = &b;
      next = &c;
      (*next)[i_seq_pivot].push_back(
        node(asu_mappings, i_seq_pivot, sgtbx::rt_mx(1,1)));
    }

    void
    shift()
    {
      container* tmp = prev;
      prev = middle;
      middle = next;
      next = tmp;
      next->clear();
      next->resize(middle->size());
    }

    unsigned
    count_next() const
    {
      unsigned result = 0;
      for(container::const_iterator i=next->begin(); i!=next->end(); i++) {
        result += i->size();
      }
      return result;
    }

    container a;
    container b;
    container c;
    container* prev;
    container* middle;
    container* next;
  };

  af::shared<std::vector<unsigned> >
  simple(
    crystal::pair_asu_table<> const& pair_asu_table,
    unsigned n_shells)
  {
    direct_space_asu::asu_mappings<> const&
      asu_mappings = *pair_asu_table.asu_mappings().get();
    af::const_ref<pair_asu_dict>
      asu_table_ref = pair_asu_table.table().const_ref();
    CCTBX_ASSERT(asu_mappings.mappings_const_ref().size()
              == asu_table_ref.size());
    af::shared<std::vector<unsigned> > term_table;
    for(unsigned i_seq_pivot=0;
                 i_seq_pivot<asu_table_ref.size();
                 i_seq_pivot++) {
      pair_asu_dict pair_asu_dict_pivot = asu_table_ref[i_seq_pivot];
      sgtbx::rt_mx rt_mx_pivot = asu_mappings.get_rt_mx(i_seq_pivot, 0);
      if (pair_asu_dict_pivot.size() == 0) {
        term_table.push_back(std::vector<unsigned>());
        continue;
      }
      three_shells shells(asu_mappings, i_seq_pivot);
      std::vector<unsigned> terms;
      terms.push_back(1);
      for(unsigned i_shell_minus_1=0;
                   i_shell_minus_1<n_shells;
                   i_shell_minus_1++) {
        shells.shift();
        for(unsigned i_seq_m=0;i_seq_m<shells.middle->size();i_seq_m++) {
          std::vector<node>& nodes_middle = (*shells.middle)[i_seq_m];
        for(unsigned i_node_m=0;i_node_m<nodes_middle.size();i_node_m++) {
          node node_m = nodes_middle[i_node_m];
          sgtbx::rt_mx rt_mx_i = asu_mappings.get_rt_mx(node_m.i_seq, 0);
          sgtbx::rt_mx rt_mx_ni = node_m.rt_mx.multiply(rt_mx_i.inverse());
          pair_asu_dict::const_iterator
            pair_asu_dict_end = asu_table_ref[node_m.i_seq].end();
          for(pair_asu_dict::const_iterator
                pair_asu_dict_i = asu_table_ref[node_m.i_seq].begin();
                pair_asu_dict_i != pair_asu_dict_end;
                pair_asu_dict_i++) {
            unsigned j_seq = pair_asu_dict_i->first;
            pair_asu_j_sym_groups const& j_sym_groups=pair_asu_dict_i->second;
            for(unsigned i_group=0; i_group<j_sym_groups.size(); i_group++) {
              pair_asu_j_sym_group j_sym_group = j_sym_groups[i_group];
              pair_asu_j_sym_group::const_iterator
                j_sym_group_end = j_sym_group.end();
              for(pair_asu_j_sym_group::const_iterator
                    j_sym_group_i = j_sym_group.begin();
                    j_sym_group_i != j_sym_group_end;
                    j_sym_group_i++) {
                unsigned j_sym = *j_sym_group_i;
                sgtbx::rt_mx rt_mx_j = asu_mappings.get_rt_mx(j_seq, j_sym);
                node new_node(asu_mappings, j_seq, rt_mx_ni.multiply(rt_mx_j));
                if (   !find_node(new_node, (*shells.prev)[j_seq])
                    && !find_node(new_node, (*shells.middle)[j_seq])
                    && !find_node(new_node, (*shells.next)[j_seq])) {
                  (*shells.next)[j_seq].push_back(new_node);
                }
              }
            }
          }
        }}
        terms.push_back(shells.count_next());
      }
      term_table.push_back(terms);
    }
    return term_table;
  }

}}} // namespace cctbx::crystal::coordination_sequences

#endif // CCTBX_CRYSTAL_COORDINATION_SEQUENCES_H
