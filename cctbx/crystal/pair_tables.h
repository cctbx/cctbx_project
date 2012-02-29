#ifndef CCTBX_CRYSTAL_PAIR_TABLES_H
#define CCTBX_CRYSTAL_PAIR_TABLES_H

#include <cctbx/crystal/neighbors_fast.h>
#include <cctbx/eltbx/covalent_radii.h>
#include <boost/scoped_array.hpp>
#include <map>
#include <set>

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
    af::const_ref<pair_sym_dict> const& pair_sym_table,
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

  class adp_iso_local_sphere_restraints_energies
  {
    public:
      adp_iso_local_sphere_restraints_energies() {}

      adp_iso_local_sphere_restraints_energies(
        af::const_ref<pair_sym_dict> const& pair_sym_table,
        scitbx::mat3<double> const& orthogonalization_matrix,
        af::const_ref<scitbx::vec3<double> > const& sites_frac,
        af::const_ref<double> const& u_isos,
        af::const_ref<bool> const& selection,
        af::const_ref<bool> const& use_u_iso,
        af::const_ref<bool> const& grad_u_iso,
        double sphere_radius,
        double distance_power,
        double average_power,
        double min_u_sum,
        bool compute_gradients,
        bool collect)
      {
        CCTBX_ASSERT(sites_frac.size() == pair_sym_table.size());
        CCTBX_ASSERT(u_isos.size() == pair_sym_table.size());
        CCTBX_ASSERT(use_u_iso.size() == pair_sym_table.size());
        CCTBX_ASSERT(grad_u_iso.size() == pair_sym_table.size());
        CCTBX_ASSERT(selection.size() == pair_sym_table.size());
        number_of_restraints = 0;
        residual_sum = 0;
        if (compute_gradients) {
          gradients.resize(u_isos.size());
        }
        for(unsigned i_seq=0;i_seq<pair_sym_table.size();i_seq++) {
          if(use_u_iso[i_seq] && selection[i_seq]) {
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
                double dist = (orthogonalization_matrix
                            * (site_i - rt_mx_list[i_op] * site_j)).length();
                if (dist <= sphere_radius && dist > 0.0) {
                  if(dist < 0.1) {
                    dist = 0.1;
                  }
                  double one_over_weight = std::pow(dist, distance_power);
                  CCTBX_ASSERT(one_over_weight != 0);
                  double weight = 1./one_over_weight;
                  double u1 = u_isos[i_seq];
                  double u2 = u_isos[j_seq];
                  if(u1 >= 0.0 && u2 >= 0) {
                     double sum = u1 + u2;
                     if (sum >= min_u_sum) {
                       double ave_pow = std::pow(sum/2, average_power);
                       CCTBX_ASSERT(ave_pow != 0);
                       double u_diff = u1 - u2;
                       double u_diff_sq = u_diff * u_diff;
                       number_of_restraints++;
                       residual_sum += (weight * u_diff_sq / ave_pow);
                       if (compute_gradients && grad_u_iso[i_seq] &&
                           grad_u_iso[j_seq]) {
                         double mem1 = 2.* u_diff / ave_pow;
                         CCTBX_ASSERT(sum * ave_pow != 0);
                         double mem2 = average_power * u_diff_sq / (sum*ave_pow);
                         gradients[i_seq] += weight * (mem1 - mem2);
                         gradients[j_seq] += weight * (-mem1 - mem2);
                       }
                       if (collect) {
                         u_i.push_back(u1);
                         u_j.push_back(u2);
                         r_ij.push_back(dist);
                     }
                   }
                }
              }
            }
          }
         }
        }
      }

      unsigned number_of_restraints;
      double residual_sum;
      af::shared<double> gradients;
      af::shared<double> u_i;
      af::shared<double> u_j;
      af::shared<double> r_ij;
  };

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
    af::const_ref<pair_sym_dict> const& pair_sym_table,
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
     typedef std::map<std::string, FloatType> RadiiRegType;
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
      {}

      //! Support for incremental_pairs. Not available from Python.
      void
      reserve(std::size_t n_sites_final) { table_.reserve(n_sites_final); }

      //! Support for incremental_pairs. Not available from Python.
      void
      grow(std::size_t number_of_additional_entries)
      {
        table_.resize(table_.size()+number_of_additional_entries);
      }

      //! Instance as passed to the constructor.
      boost::shared_ptr<
        direct_space_asu::asu_mappings<FloatType, IntShiftType> >
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

      //! Pair count for each site.
      af::shared<std::size_t>
      pair_counts() const
      {
        af::const_ref<pair_asu_dict> tab = table_.const_ref();
        af::shared<std::size_t> result((af::reserve(tab.size())));
        for(unsigned i_seq=0;i_seq<tab.size();i_seq++) {
          std::size_t count = 0;
          pair_asu_dict const& dict = tab[i_seq];
          for(pair_asu_dict::const_iterator
                dict_i=dict.begin();
                dict_i!=dict.end();
                dict_i++) {
            pair_asu_j_sym_groups const& jgs = dict_i->second;
            for(unsigned ig=0;ig<jgs.size();ig++) {
              count += jgs[ig].size();
            }
          }
          result.push_back(count);
        }
        return result;
      }

      //! Returns selection for cluster pivots.
      /*! The list of sites is assumed to be ordered by significance
          (e.g. peak height), most significant first.
       */
      af::shared<std::size_t>
      cluster_pivot_selection(
        bool general_positions_only=false,
        std::size_t max_clusters=0,
        unsigned estimated_reduction_factor=4) const
      {
        af::const_ref<pair_asu_dict> tab = table_.const_ref();
        af::shared<std::size_t> result;
        if (max_clusters > 0) {
          result.reserve(max_clusters);
        }
        else if (estimated_reduction_factor > 1) {
          result.reserve(  (tab.size()+estimated_reduction_factor-1)
                         / estimated_reduction_factor);
        }
        else {
          result.reserve(tab.size());
        }
        boost::scoped_array<bool> keep_flags(new bool[tab.size()]);
        for(unsigned i_seq=0;i_seq<tab.size();i_seq++) {
          if (general_positions_only
              && asu_mappings_->site_symmetry_table()
                   .is_special_position(i_seq)) {
            keep_flags[i_seq] = false;
            continue;
          }
          bool keep = true;
          pair_asu_dict const& dict = tab[i_seq];
          for(pair_asu_dict::const_iterator
                dict_i=dict.begin();
                dict_i!=dict.end();
                dict_i++) {
            unsigned j_seq = dict_i->first;
            if (j_seq == i_seq || (j_seq < i_seq && keep_flags[j_seq])) {
              keep = false;
              break;
            }
          }
          keep_flags[i_seq] = keep;
          if (keep) {
            result.push_back(i_seq);
            if (result.size() == max_clusters) break;
          }
        }
        return result;
      }

      /*! \brief Uses neighbors::fast_pair_generator to add all pairs with
          distances <= distance_cutoff*(1+epsilon).
       */
      /*! All symmetrically equivalent pairs are automatically generated.
       */
      pair_asu_table&
      add_all_pairs(
        FloatType const& distance_cutoff,
        FloatType const& min_cubicle_edge=5,
        FloatType const& epsilon=1.e-6)
      {
        neighbors::fast_pair_generator<FloatType, IntShiftType> pair_generator(
          asu_mappings_owner_,
          distance_cutoff*(1+epsilon),
          /*minimal*/ true,
          min_cubicle_edge);
        while (!pair_generator.at_end()) {
          add_pair(pair_generator.next());
        }
        return *this;
      }

      /*! \brief Uses neighbors::fast_pair_generator to add all pairs with
          distances <= [sum(covalent radii) + tolerance].
       */
      /*! All symmetrically equivalent pairs are automatically generated.
       */
      pair_asu_table&
      add_covalent_pairs(
        af::const_ref<std::string> const& scattering_types,
        af::const_ref<std::string> const& exclude_scattering_types
          = af::const_ref<std::string>(0,0),
        af::const_ref<std::size_t> const& conformer_indices
          = af::const_ref<std::size_t>(0, 0),
        af::const_ref<std::size_t> const& sym_excl_indices
          = af::const_ref<std::size_t>(0, 0),
        FloatType const& distance_cutoff=3.5,
        FloatType const& min_cubicle_edge=5,
        FloatType const& tolerance=0.5,
        FloatType const& epsilon=1.e-6,
        RadiiRegType const &radii=RadiiRegType())
      {
        CCTBX_ASSERT(!conformer_indices.size()
                  ||  conformer_indices.size() == scattering_types.size());
        CCTBX_ASSERT(!sym_excl_indices.size()
                  ||  sym_excl_indices.size() == scattering_types.size());
        neighbors::fast_pair_generator<FloatType, IntShiftType> pair_generator(
          asu_mappings_owner_,
          distance_cutoff*(1+epsilon),
          /*minimal*/ true,
          min_cubicle_edge);
        RadiiRegType radii_;
        for (std::size_t i=0; i < scattering_types.size(); i++) {
          radii_[scattering_types[i]] =
            eltbx::covalent_radii::table(scattering_types[i]).radius();
        }
        typedef typename RadiiRegType::const_iterator Itr;
        for (Itr itr = radii.begin(); itr != radii.end(); itr++) {
          radii_[itr->first] = itr->second;
        }

        while (!pair_generator.at_end()) {
          direct_space_asu::asu_mapping_index_pair_and_diff<FloatType>
            const& pair = pair_generator.next();
          if (   std::find(exclude_scattering_types.begin(),
                           exclude_scattering_types.end(),
                           scattering_types[pair.i_seq])
                   != exclude_scattering_types.end()
              || std::find(exclude_scattering_types.begin(),
                           exclude_scattering_types.end(),
                           scattering_types[pair.j_seq])
                   != exclude_scattering_types.end()) {
               continue;
          }
          if (   conformer_indices.size()
              && conformer_indices[pair.i_seq] != 0
              && conformer_indices[pair.j_seq] != 0
              && conformer_indices[pair.i_seq]
              != conformer_indices[pair.j_seq]) {
                continue;
          }
          if (   sym_excl_indices.size()
              && sym_excl_indices[pair.i_seq] != 0
              && sym_excl_indices[pair.j_seq] != 0
              && asu_mappings_->get_rt_mx_i(pair)
              != asu_mappings_->get_rt_mx_j(pair)) {
                continue; // don't bond to sym equivs
          }
          if (   conformer_indices.size()
              && sym_excl_indices.size()
              && ((   conformer_indices[pair.i_seq] != 0
                   && sym_excl_indices[pair.j_seq] != 0)
                 ||
                  (   conformer_indices[pair.j_seq] != 0
                   && sym_excl_indices[pair.i_seq] != 0))) {
                continue;
          }
          FloatType const max_bond_length = radii_[scattering_types[pair.i_seq]] +
            radii_[scattering_types[pair.j_seq]] + tolerance;
          if (std::sqrt(pair.dist_sq) <= max_bond_length) {
            add_pair(pair);
          }
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
        return *this;
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

      //! Extracts pair_sym_table of interaction pairs.
      /*! The default is to consider only pairs unique under symmetry.
          The result may then be interpreted as an "asymmetric unit of pair
          interactions." The returned pairs read (i_seq, j_seq, rt_mx_ji)
          where rt_mx_ji is the symmetry moving j_seq to the site interacting
          with i_seq. Pairs such that j_seq is in the asu, i.e. rt_mx_ji is the
          unit Seitz matrix are preferred.

          If skip_j_seq_less_than_i_seq == false the result may
          contain redundancies. This option is useful for debugging
          and development purposes but also when one needs for each site to
          have directly access to all the sites it interacts with.

          If all_interactions_from_inside_asu is true, all the interactions
          involving any site in the asu are generated, even if some of them are
          symmetry-equivalent. The typical case is that of a site C1
          on a special position which is bonded to a site C2' that is the
          symmetry equivalent of a site C2 in the asu through an operator
          leaving C1 invariant. Then in addition to the interaction C1--C2',
          there is the interaction C1--C2. When all_interactions_from_inside_asu
          is false, the latter interaction will be generated. When that flag
          is true, the both of C1--C2' and C1--C2 will be generated.
       */
      pair_sym_table
      extract_pair_sym_table(bool skip_j_seq_less_than_i_seq=true,
                             bool all_interactions_from_inside_asu=false) const
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
            sgtbx::site_symmetry_ops const
            &site_symm_j_seq = asu_mappings()->site_symmetry_table().get(j_seq);
            std::vector<sgtbx::rt_mx>& rt_mx_list = sym_dict[j_seq];
            for(pair_asu_j_sym_groups::const_iterator
                  j_sym_groups_i = asu_dict_i->second.begin();
                  j_sym_groups_i != asu_dict_i->second.end();
                  j_sym_groups_i++) {
              boost::optional<sgtbx::rt_mx> first_rt_mx_ji, candidate_rt_mx_ji;
              for (pair_asu_j_sym_group::const_iterator
                   j_sym_i = j_sym_groups_i->begin();
                   j_sym_i != j_sym_groups_i->end();
                   ++j_sym_i) {
                unsigned j_sym = *j_sym_i;
                sgtbx::rt_mx rt_mx_j = asu_mappings_->get_rt_mx(j_seq, j_sym);
                sgtbx::rt_mx rt_mx_ji = rt_mx_i_inv.multiply(rt_mx_j);
                if (!first_rt_mx_ji) first_rt_mx_ji = rt_mx_ji;
                if (site_symm_j_seq.contains(rt_mx_ji)) {
                  rt_mx_ji = rt_mx_ji.unit_mx();
                  if (!all_interactions_from_inside_asu) {
                    candidate_rt_mx_ji = rt_mx_ji;
                    break;
                  }
                }
                if (all_interactions_from_inside_asu) {
                  rt_mx_list.push_back(rt_mx_ji);
                }
              }
              if (all_interactions_from_inside_asu) continue;
              if (candidate_rt_mx_ji) rt_mx_list.push_back(*candidate_rt_mx_ji);
              else if (first_rt_mx_ji) rt_mx_list.push_back(*first_rt_mx_ji);
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
          boost::optional<sgtbx::rt_mx> rt_mx_ji_inv;
          if (i_seq == j_seq) {
            rt_mx_ji_inv = rt_mx_ji.inverse();
          }
          for(unsigned i_mi=0;i_mi<site_symmetry_matrices.size();i_mi++) {
            sgtbx::rt_mx const& mi = site_symmetry_matrices[i_mi];
            if (rt_mx_ji_inv) {
              sgtbx::rt_mx rt_mx_j_eq
                = rt_mx_i.multiply(mi.multiply(*rt_mx_ji_inv));
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
      }

      //! Constructs new pair_asu_table with interactions for all bond angles.
      pair_asu_table
      angle_pair_asu_table() const
      {
        pair_asu_table result(asu_mappings_owner_);
        af::const_ref<pair_asu_dict> table_ref = table_.const_ref();
        typedef af::tiny<unsigned, 2> tu2;
        std::vector<tu2> pair_list;
        pair_list.reserve(16); // very generous estimate
        for(unsigned i_seq=0;i_seq<table_.size();i_seq++) {
          pair_list.clear();
          pair_asu_dict const& asu_dict = table_ref[i_seq];
          for(pair_asu_dict::const_iterator
                asu_dict_i = asu_dict.begin();
                asu_dict_i != asu_dict.end();
                asu_dict_i++) {
            unsigned j_seq = asu_dict_i->first;
            for(pair_asu_j_sym_groups::const_iterator
                  j_sym_groups_i = asu_dict_i->second.begin();
                  j_sym_groups_i != asu_dict_i->second.end();
                  j_sym_groups_i++) {
              for(pair_asu_j_sym_group::const_iterator
                    j_sym_i = j_sym_groups_i->begin();
                    j_sym_i != j_sym_groups_i->end();
                    j_sym_i++) {
                pair_list.push_back(tu2(j_seq, *j_sym_i));
              }
            }
          }
          unsigned pair_list_size = static_cast<unsigned>(pair_list.size());
          for(unsigned i_jj1=0;i_jj1+1<pair_list_size;i_jj1++) {
            tu2 const& jj1 = pair_list[i_jj1];
            sgtbx::rt_mx rt_mx_jj1_inv = asu_mappings_->get_rt_mx(
              jj1[0], jj1[1]).inverse();
            for(unsigned i_jj2=i_jj1+1;i_jj2<pair_list_size;i_jj2++) {
              tu2 const& jj2 = pair_list[i_jj2];
              result.add_pair(
                jj1[0],
                jj2[0],
                rt_mx_jj1_inv.multiply(asu_mappings_->get_rt_mx(
                  jj2[0], jj2[1])));
            }
          }
        }
        return result;
      }

    protected:
      boost::shared_ptr<
        direct_space_asu::asu_mappings<FloatType, IntShiftType> >
          asu_mappings_owner_;
      const direct_space_asu::asu_mappings<FloatType, IntShiftType>*
        asu_mappings_;
      pair_asu_table_table table_;
  };

}} // namespace cctbx::crystal

#endif // CCTBX_CRYSTAL_PAIR_TABLES_H
