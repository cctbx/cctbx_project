#ifndef CCTBX_CRYSTAL_PAIR_TABLES_H
#define CCTBX_CRYSTAL_PAIR_TABLES_H

#include <cctbx/crystal/neighbors_fast.h>
#include <map>

namespace cctbx { namespace crystal {

  typedef sgtbx::rt_mx pair_sym_op;
  typedef std::vector<sgtbx::rt_mx> pair_sym_ops;
  typedef std::map<unsigned, pair_sym_ops> pair_sym_dict;
  typedef af::shared<pair_sym_dict> pair_sym_table;

  typedef std::set<unsigned> pair_asu_j_sym_group;
  typedef std::vector<pair_asu_j_sym_group> pair_asu_j_sym_groups;
  typedef std::map<unsigned, pair_asu_j_sym_groups> pair_asu_dict;
  typedef af::shared<pair_asu_dict> pair_asu_table_table;

  template <typename FloatType=double, typename IntShiftType=int>
  class pair_asu_table
  {
    public:
      pair_asu_table() {}

      pair_asu_table(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<
            FloatType, IntShiftType> >& asu_mappings)
      :
        asu_mappings_owner_(asu_mappings),
        asu_mappings_(asu_mappings.get()),
        table_(asu_mappings->mappings_const_ref().size())
      {
        asu_mappings->lock();
      }

      //! Instance as passed to the constructor.
      boost::shared_ptr<direct_space_asu::asu_mappings<> > const&
      asu_mappings() const { return asu_mappings_owner_; }

      pair_asu_table_table const&
      table() const { return table_; }

      bool
      contains(direct_space_asu::asu_mapping_index_pair const& pair)
      {
        return contains(pair.i_seq, pair.j_seq, pair.j_sym);
      }

      bool
      contains(unsigned i_seq, unsigned j_seq, unsigned j_sym)
      {
        pair_asu_dict const& asu_dict = table_[i_seq];
        pair_asu_dict::const_iterator asu_dict_i = asu_dict.find(j_seq);
        if (asu_dict_i != asu_dict.end()) {
          for (pair_asu_j_sym_groups::const_iterator
                 j_sym_groups_i = asu_dict_i->second.begin();
                 j_sym_groups_i != asu_dict_i->second.end();
                 j_sym_groups_i++) {
            if (j_sym_groups_i->find(j_sym) != j_sym_groups_i->end()) {
              return true;
            }
          }
        }
        return false;
      }

      pair_asu_table&
      add_all_pairs(
        FloatType const& distance_cutoff,
        FloatType const& epsilon=1.e-6)
      {
        bool minimal = true;
        neighbors::fast_pair_generator<FloatType, IntShiftType> pair_generator(
          asu_mappings_owner_,
          distance_cutoff,
          minimal);
        while (!pair_generator.at_end()) {
          direct_space_asu::asu_mapping_index_pair
            pair = pair_generator.next();
          sgtbx::rt_mx rt_mx_i = asu_mappings_->get_rt_mx_i(pair);
          sgtbx::rt_mx rt_mx_j = asu_mappings_->get_rt_mx_j(pair);
          add_pair(
            pair.i_seq,
            pair.j_seq,
            rt_mx_i.inverse().multiply(rt_mx_j));
        }
        return *this;
      }

      pair_asu_table&
      add_pair_sym_table(pair_sym_table const& sym_table)
      {
        af::const_ref<pair_sym_dict> sym_table_ref = sym_table.const_ref();
        for(unsigned i_seq=0;i_seq<sym_table_ref.size();i_seq++) {
          for(pair_sym_dict::const_iterator
                sym_dict_i = sym_table_ref[i_seq].begin();
                sym_dict_i != sym_table_ref[i_seq].end();
                sym_dict_i++) {
            unsigned j_seq = sym_dict_i->first;
            std::vector<sgtbx::rt_mx> const& rt_mx_list = sym_dict_i->second;
            for(unsigned i=0;i<rt_mx_list.size();i++) {
              add_pair(i_seq, j_seq, rt_mx_list[i]);
            }
          }
        }
        return *this;
      }

      pair_asu_table&
      add_pair(unsigned i_seq, unsigned j_seq, sgtbx::rt_mx const& rt_mx_ji)
      {
        bool is_new = process_pair(i_seq, j_seq, rt_mx_ji);
        if (is_new && i_seq != j_seq) {
          is_new = process_pair(j_seq, i_seq, rt_mx_ji.inverse_cancel());
          CCTBX_ASSERT(is_new);
        }
        return *this;
      }

      pair_asu_table&
      add_pair(af::tiny<unsigned, 2> const& i_seqs)
      {
        sgtbx::rt_mx rt_mx_ji(1, 1);
        bool is_new = process_pair(i_seqs[0], i_seqs[1], rt_mx_ji);
        if (is_new && i_seqs[0] != i_seqs[1]) {
          is_new = process_pair(i_seqs[1], i_seqs[0], rt_mx_ji);
          CCTBX_ASSERT(is_new);
        }
        return *this;
      }

      pair_sym_table
      extract_pair_sym_table() const
      {
        pair_sym_table sym_table(asu_mappings_->mappings_const_ref().size());
        af::const_ref<pair_asu_dict> table_ref = table_.const_ref();
        for(unsigned i_seq=0;i_seq<table_.size();i_seq++) {
          sgtbx::rt_mx
            rt_mx_i_inv = asu_mappings_->get_rt_mx(i_seq, 0).inverse();
          pair_asu_dict const& asu_dict = table_ref[i_seq];
          pair_sym_dict& sym_dict = sym_table[i_seq];
          for(pair_asu_dict::const_iterator
                asu_dict_i = asu_dict.begin();
                asu_dict_i != asu_dict.end();
                asu_dict_i++) {
            unsigned j_seq = asu_dict_i->first;
            if (j_seq < i_seq) continue;
            std::vector<sgtbx::rt_mx>& rt_mx_list = sym_dict[j_seq];
            for(pair_asu_j_sym_groups::const_iterator
                  j_sym_groups_i = asu_dict_i->second.begin();
                  j_sym_groups_i != asu_dict_i->second.end();
                  j_sym_groups_i++) {
              unsigned j_sym = *(j_sym_groups_i->begin());
              sgtbx::rt_mx rt_mx_j = asu_mappings_->get_rt_mx(j_seq, j_sym);
              rt_mx_list.push_back(rt_mx_i_inv.multiply(rt_mx_j));
            }
          }
        }
        return sym_table;
      }

    protected:
      bool
      process_pair(
        unsigned i_seq,
        unsigned j_seq,
        sgtbx::rt_mx const& rt_mx_ji)
      {
        sgtbx::rt_mx rt_mx_i = asu_mappings_->get_rt_mx(i_seq, 0);
        sgtbx::rt_mx rt_mx_j = rt_mx_i.multiply(rt_mx_ji);
        int j_sym = asu_mappings_->find_i_sym(j_seq, rt_mx_j);
        CCTBX_ASSERT(j_sym >= 0);
        if (contains(i_seq, j_seq, j_sym)) {
          return false;
        }
        table_[i_seq][j_seq].push_back(pair_asu_j_sym_group());
        pair_asu_j_sym_group& j_syms = table_[i_seq][j_seq].back();
        af::const_ref<sgtbx::rt_mx> const&
          site_symmetry_matrices = asu_mappings_->site_symmetry_table()
            .get(i_seq).matrices().const_ref();
        for(unsigned i_mi=0;i_mi<site_symmetry_matrices.size();i_mi++) {
          sgtbx::rt_mx const& mi = site_symmetry_matrices[i_mi];
          if (i_seq == j_seq) {
            sgtbx::rt_mx rt_mx_j_eq
              = rt_mx_i.multiply(rt_mx_ji.multiply(mi).inverse_cancel());
            int j_sym_eq = asu_mappings_->find_i_sym(j_seq, rt_mx_j_eq);
            CCTBX_ASSERT(j_sym_eq >= 0);
            j_syms.insert(j_sym_eq);
          }
          sgtbx::rt_mx rt_mx_j_eq = rt_mx_i.multiply(mi.multiply(rt_mx_ji));
          int j_sym_eq = asu_mappings_->find_i_sym(j_seq, rt_mx_j_eq);
          CCTBX_ASSERT(j_sym_eq >= 0);
          j_syms.insert(j_sym_eq);
        }
        return true;
      };

      boost::shared_ptr<direct_space_asu::asu_mappings<> > asu_mappings_owner_;
      const direct_space_asu::asu_mappings<>* asu_mappings_;
      pair_asu_table_table table_;
  };

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_PAIR_TABLES_H
