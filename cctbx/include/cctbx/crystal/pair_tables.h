#ifndef CCTBX_CRYSTAL_PAIR_TABLES_H
#define CCTBX_CRYSTAL_PAIR_TABLES_H

#include <cctbx/crystal/neighbors_fast.h>
#include <map>

namespace cctbx { namespace crystal {

  //! Symmetry operation characterizing one pair interaction.
  typedef sgtbx::rt_mx pair_sym_op;
  //! Group of symmetry operations for a given i_seq and j_seq.
  typedef std::vector<sgtbx::rt_mx> pair_sym_ops;
  //! Dictionary of pair interactions for a given i_seq.
  typedef std::map<unsigned, pair_sym_ops> pair_sym_dict;
  //! Table of pair interactions indexed by i_seq.
  typedef af::shared<pair_sym_dict> pair_sym_table;

  /*! \brief Determination of distances of all pair interactions
      defined by pair_sym_table.
   */
  inline
  af::shared<double>
  get_distances(
    af::const_ref<crystal::pair_sym_dict> const& pair_sym_table,
    scitbx::mat3<double> const& orthogonalization_matrix,
    af::const_ref<scitbx::vec3<double> > const& sites_frac)
  {
    CCTBX_ASSERT(sites_frac.size() == pair_sym_table.size());
    af::shared<double> distances;
    for(unsigned i_seq=0;i_seq<pair_sym_table.size();i_seq++) {
      crystal::pair_sym_dict const& pair_sym_dict = pair_sym_table[i_seq];
      scitbx::vec3<double> const& site_i = sites_frac[i_seq];
      for(crystal::pair_sym_dict::const_iterator
            pair_sym_dict_i=pair_sym_dict.begin();
            pair_sym_dict_i!=pair_sym_dict.end();
            pair_sym_dict_i++) {
        unsigned j_seq = pair_sym_dict_i->first;
        scitbx::vec3<double> const& site_j = sites_frac[j_seq];
        af::const_ref<sgtbx::rt_mx> rt_mx_list = af::make_const_ref(
          pair_sym_dict_i->second);
        for(unsigned i_op=0;i_op<rt_mx_list.size();i_op++) {
          distances.push_back((
            orthogonalization_matrix
              * (site_i - rt_mx_list[i_op] * site_j)).length());
        }
      }
    }
    return distances;
  }

  /*! \brief Determination of distances of all pair interactions
      defined by pair_sym_table for exclusively non-crystallographic
      interactions.
   */
  /*! All symmetry operations must be the identity. An exception
      is thrown otherwise.
   */
  inline
  af::shared<double>
  get_distances(
    af::const_ref<crystal::pair_sym_dict> const& pair_sym_table,
    af::const_ref<scitbx::vec3<double> > const& sites_cart)
  {
    CCTBX_ASSERT(sites_cart.size() == pair_sym_table.size());
    af::shared<double> distances;
    for(unsigned i_seq=0;i_seq<pair_sym_table.size();i_seq++) {
      crystal::pair_sym_dict const& pair_sym_dict = pair_sym_table[i_seq];
      scitbx::vec3<double> const& site_i = sites_cart[i_seq];
      for(crystal::pair_sym_dict::const_iterator
            pair_sym_dict_i=pair_sym_dict.begin();
            pair_sym_dict_i!=pair_sym_dict.end();
            pair_sym_dict_i++) {
        unsigned j_seq = pair_sym_dict_i->first;
        scitbx::vec3<double> const& site_j = sites_cart[j_seq];
        af::const_ref<sgtbx::rt_mx> rt_mx_list = af::make_const_ref(
          pair_sym_dict_i->second);
        CCTBX_ASSERT(rt_mx_list.size() == 1);
        CCTBX_ASSERT(rt_mx_list[0].is_unit_mx());
        distances.push_back((site_i - site_j).length());
      }
    }
    return distances;
  }

  //! Set of j_sym indices of symmetrically equivalent pair interactions.
  typedef std::set<unsigned> pair_asu_j_sym_group;
  //! Array of sets of symmetrically equivalent pair interactions.
  typedef std::vector<pair_asu_j_sym_group> pair_asu_j_sym_groups;
  //! Dictionary of pair interactions for a given i_seq.
  typedef std::map<unsigned, pair_asu_j_sym_groups> pair_asu_dict;
  //! Table of pair interactions indexed by i_seq.
  typedef af::shared<pair_asu_dict> pair_asu_table_table;

  /*! \brief Managed table of pair interactions based on
      direct_space_asu::asu_mappings.
   */
  template <typename FloatType=double, typename IntShiftType=int>
  class pair_asu_table
  {
    public:
      //! Default constructor. Some data members are not initialized!
      pair_asu_table() {}

      //! Initialization of an empty table().
      pair_asu_table(
        boost::shared_ptr<
          direct_space_asu::asu_mappings<
            FloatType, IntShiftType> > asu_mappings)
      :
        asu_mappings_owner_(asu_mappings),
        asu_mappings_(asu_mappings.get()),
        table_(asu_mappings->mappings_const_ref().size())
      {
        asu_mappings->lock();
      }

      //! Instance as passed to the constructor.
      boost::shared_ptr<direct_space_asu::asu_mappings<> >
      asu_mappings() const { return asu_mappings_owner_; }

      //! Access to raw table.
      pair_asu_table_table const&
      table() const { return table_; }

      //! True if pair is in table().
      bool
      contains(direct_space_asu::asu_mapping_index_pair const& pair) const
      {
        return contains(pair.i_seq, pair.j_seq, pair.j_sym);
      }

      //! True if pair characterized by i_seq, j_seq, j_sym is in table().
      bool
      contains(unsigned i_seq, unsigned j_seq, unsigned j_sym) const
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

      //! The order of the pair_asu_j_sym_groups is irrelevant.
      bool
      operator==(pair_asu_table const& other) const
      {
        af::const_ref<pair_asu_dict> tab_a = table_.const_ref();
        af::const_ref<pair_asu_dict> tab_b = other.table_.const_ref();
        if (tab_a.size() != tab_b.size()) return false;
        for(unsigned i_seq=0;i_seq<tab_a.size();i_seq++) {
          pair_asu_dict const& dict_a = tab_a[i_seq];
          pair_asu_dict const& dict_b = tab_b[i_seq];
          if (dict_a.size() != dict_b.size()) return false;
          for(pair_asu_dict::const_iterator
                dict_a_i=dict_a.begin();
                dict_a_i!=dict_a.end();
                dict_a_i++) {
            pair_asu_dict::const_iterator
              dict_b_i = dict_b.find(dict_a_i->first);
            if (dict_b_i == dict_b.end()) return false;
            pair_asu_j_sym_groups const& jgs_a = dict_a_i->second;
            pair_asu_j_sym_groups const& jgs_b = dict_b_i->second;
            if (jgs_a.size() != jgs_b.size()) return false;
            std::vector<bool> is_matched(jgs_b.size(), false);
            for(unsigned ig_a=0;ig_a<jgs_a.size();ig_a++) {
              unsigned ig_b=0;
              for(;ig_b<jgs_b.size();ig_b++) {
                if (is_matched[ig_b]) continue;
                if (jgs_a[ig_a] == jgs_b[ig_b]) {
                  is_matched[ig_b] = true;
                  break;
                }
              }
              if (ig_b == jgs_b.size()) return false;
            }
          }
        }
        return true;
      }

      //! Shorthand for: !operator==
      bool
      operator!=(pair_asu_table const& other) const
      {
        return !((*this) == other);
      }

      /*! \brief Uses neighbors::fast_pair_generator to add all pairs with
          distances <= distance_cutoff*(1+epsilon).
       */
      /*! All symmetrically equivalent pairs are automatically generated.
       */
      pair_asu_table&
      add_all_pairs(
        FloatType const& distance_cutoff,
        FloatType const& epsilon=1.e-6)
      {
        bool minimal = true;
        neighbors::fast_pair_generator<FloatType, IntShiftType> pair_generator(
          asu_mappings_owner_,
          distance_cutoff*(1+epsilon),
          minimal);
        while (!pair_generator.at_end()) {
          add_pair(pair_generator.next());
        }
        return *this;
      }

      //! Adds all pairs defined by sym_table.
      /*! All symmetrically equivalent pairs are automatically generated.
       */
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

      //! Adds a pair and all symmetrically equivalent pairs.
      pair_asu_table&
      add_pair(direct_space_asu::asu_mapping_index_pair const& pair)
      {
        sgtbx::rt_mx rt_mx_i = asu_mappings_->get_rt_mx_i(pair);
        sgtbx::rt_mx rt_mx_j = asu_mappings_->get_rt_mx_j(pair);
        add_pair(
          pair.i_seq,
          pair.j_seq,
          rt_mx_i.inverse().multiply(rt_mx_j));
      }

      //! Adds the pair defined by i_seq, j_seq and rt_mx_ji.
      /*! rt_mx_ji is the symmetry operation that maps the original site
          referenced by j_seq to the site interacting with the original
          site referenced by i_seq.

          All symmetrically equivalent pairs are automatically generated.
       */
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

      /*! \brief Adds the pair defined by i_seq, j_seq assuming
          rt_mx_ji is the identity matrix.
       */
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

      //! Extracts pair_sym_table of interactions unique under symmetry.
      /*! The result may be interpreted as an "asymmetric unit of pair
          interactions."

          If skip_j_seq_less_than_i_seq == false the result may
          contain redundancies. This option is mainly for debugging
          and development purposes.
       */
      pair_sym_table
      extract_pair_sym_table(bool skip_j_seq_less_than_i_seq=true) const
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
            if (skip_j_seq_less_than_i_seq && j_seq < i_seq) continue;
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

      /*! \brief Addition of a pair interaction and all its
          symmetrically equivalent interactions to table().
       */
      /*! This functions calls the other overload to perform
          the actual work.

          Not available in Python.
       */
      bool
      process_pair(
        unsigned i_seq,
        unsigned j_seq,
        sgtbx::rt_mx const& rt_mx_ji)
      {
        sgtbx::rt_mx rt_mx_i = asu_mappings_->get_rt_mx(i_seq, 0);
        sgtbx::rt_mx rt_mx_j = rt_mx_i.multiply(rt_mx_ji);
        int j_sym = asu_mappings_->find_i_sym(j_seq, rt_mx_j);
        return process_pair(i_seq, j_seq, rt_mx_ji, rt_mx_i, j_sym);
      }

      /*! \brief Addition of a pair interaction and all its
          symmetrically equivalent interactions to table().
       */
      /*! The site_symmetry of i_seq is used to generate the
          symmetrically equivalent interactions.

          Not available in Python.
       */
      bool
      process_pair(
        unsigned i_seq,
        unsigned j_seq,
        sgtbx::rt_mx const& rt_mx_ji,
        sgtbx::rt_mx const& rt_mx_i,
        int j_sym)
      {
        CCTBX_ASSERT(j_sym >= 0);
        if (contains(i_seq, j_seq, j_sym)) {
          return false;
        }
        table_[i_seq][j_seq].push_back(pair_asu_j_sym_group());
        pair_asu_j_sym_group& j_syms = table_[i_seq][j_seq].back();
        sgtbx::site_symmetry_table const&
          site_symmetry_table = asu_mappings_->site_symmetry_table();
        if (   i_seq != j_seq
            && !site_symmetry_table.is_special_position(i_seq)) {
          j_syms.insert(j_sym);
        }
        else {
          af::const_ref<sgtbx::rt_mx> const&
            site_symmetry_matrices = site_symmetry_table
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
        }
        return true;
      };

    protected:
      boost::shared_ptr<direct_space_asu::asu_mappings<> > asu_mappings_owner_;
      const direct_space_asu::asu_mappings<>* asu_mappings_;
      pair_asu_table_table table_;
  };

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_PAIR_TABLES_H
