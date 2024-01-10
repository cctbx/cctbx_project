#include <iotbx/pdb/hierarchy.h>
#include <iotbx/pdb/common_residue_names.h>
#include <iotbx/error.h>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb { namespace hierarchy {

namespace detail {

  struct cmp_atom_labels_functor
  {
    std::vector<std::string> const& labels;

    cmp_atom_labels_functor(
      std::vector<std::string> const& labels_)
    :
      labels(labels_)
    {}

    bool
    operator()(
      unsigned i,
      unsigned j) const
    {
      return (labels[i] < labels[j]);
    }
  };

  bool
  cmp_first_element_of_vectors(
    std::vector<unsigned> const& a,
    std::vector<unsigned> const& b)
  {
    return (a[0] < b[0]);
  }

  unsigned
  find_duplicate_atom_labels(
    af::shared<af::shared<atom> >& duplicate_atom_labels,
    hierarchy::model const& model,
    unsigned model_atoms_size,
    std::vector<std::string> const& model_atom_labels)
  {
    if (model_atoms_size == 0) return 0;
    boost::scoped_array<unsigned> indices(new unsigned[model_atoms_size]);
    for(unsigned j=0;j<model_atoms_size;j++) {
      indices[j] = j;
    }
    std::sort(
      indices.get(),
      indices.get()+model_atoms_size,
      cmp_atom_labels_functor(model_atom_labels));
    std::vector<std::vector<unsigned> > groups;
    std::vector<unsigned> group;
    unsigned j_start = 0;
    for(unsigned j=1;j<model_atoms_size+1U;j++) {
      if (j != model_atoms_size) {
        if (model_atom_labels[indices[j_start]] == model_atom_labels[indices[j]]) continue;
      }
      if (j_start+1U == j) {
        j_start++;
      }
      else {
        group.reserve(j-j_start);
        while (j_start != j) group.push_back(indices[j_start++]);
        std::sort(group.begin(), group.end());
        groups.push_back(std::vector<unsigned>());
        groups.back().swap(group);
      }
    }
    unsigned result = 0;
    if (groups.size() != 0) {
      std::sort(groups.begin(), groups.end(), cmp_first_element_of_vectors);
      af::shared<atom> atoms_owner = model.atoms();
      af::const_ref<atom> atoms = atoms_owner.const_ref();
      IOTBX_ASSERT(atoms.size() == model_atoms_size);
      unsigned n_groups = static_cast<unsigned>(groups.size());
      duplicate_atom_labels.reserve(duplicate_atom_labels.size() + n_groups);
      for(unsigned i_g=0;i_g<n_groups;i_g++) {
        std::vector<unsigned> const& g = groups[i_g];
        af::shared<atom> ga((af::reserve(g.size())));
        unsigned n_i = static_cast<unsigned>(g.size());
        for(unsigned i_i=0;i_i<n_i;i_i++) {
          ga.push_back(atoms[g[i_i]]);
        }
        duplicate_atom_labels.push_back(ga);
        result += n_i;
      }
    }
    return result;
  }

} // namespace detail

  overall_counts::overall_counts(
    hierarchy::root const& root_)
  :
    root(root_),
    n_empty_models(0),
    n_empty_chains(0),
    n_empty_residue_groups(0),
    n_empty_atom_groups(0),
    n_duplicate_model_ids(0),
    n_duplicate_chain_ids(0),
    n_duplicate_atom_labels(0),
    n_models(root.models_size()),
    n_chains(0),
    n_alt_conf(0),
    n_residues(0),
    n_residue_groups(0),
    n_explicit_chain_breaks(0),
    n_atoms(0),
    n_anisou(0),
    n_alt_conf_none(0),
    n_alt_conf_pure(0),
    n_alt_conf_proper(0),
    n_alt_conf_improper(0),
    n_chains_with_mix_of_proper_and_improper_alt_conf(0)
  {
    for(unsigned i_md=0;i_md<n_models;i_md++) {
      hierarchy::model const& model = root.models()[i_md];
      if (model.chains_size() == 0) n_empty_models++;
      model_ids[model.data->id]++;
      std::map<std::string, unsigned> model_chain_ids;
      unsigned model_atoms_size = model.atoms_size();
      std::vector<std::string > model_atom_labels;
      unsigned i_model_atom = 0;
      unsigned n_ch = model.chains_size();
      n_chains += n_ch;
      for(unsigned i_ch=0;i_ch<n_ch;i_ch++) {
        hierarchy::chain const& chain = model.chains()[i_ch];
        if (chain.residue_groups_size() == 0) n_empty_chains += 1;
        model_chain_ids[chain.data->id]++;
        chain_ids[chain.data->id]++;
        std::set<std::string> chain_altlocs;
        boost::optional<residue_group> chain_alt_conf_proper;
        boost::optional<residue_group> chain_alt_conf_improper;
        bool suppress_chain_break = true;
        boost::optional<residue_group> prev_rg;
        unsigned n_rg = chain.residue_groups_size();
        for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
          residue_group const& rg = chain.residue_groups()[i_rg];
          if (rg.atom_groups_size() == 0) n_empty_residue_groups++;
          if (!rg.data->link_to_previous && !suppress_chain_break) {
            n_explicit_chain_breaks++;
          }
          suppress_chain_break = false;
          bool have_main_conf = false;
          bool have_blank_altloc = false;
          std::set<std::string> rg_altlocs;
          std::set<std::string> rg_resnames;
          std::map<std::string, std::vector<std::string> > altloc_resnames;
          unsigned n_ag = rg.atom_groups_size();
          for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
            atom_group const& ag = rg.atom_groups()[i_ag];
            if (ag.atoms_size() == 0) n_empty_atom_groups++;
            std::string altloc = ag.data->altloc;
            if (altloc == "") {
              have_main_conf = true;
            }
            else {
              if (altloc == blank_altloc_string) {
                have_blank_altloc = true;
              }
              rg_altlocs.insert(altloc);
            }
            rg_resnames.insert(ag.data->resname);
            altloc_resnames[altloc].push_back(ag.data->resname);
            unsigned n_ats = ag.atoms_size();
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              hierarchy::atom const& atom = ag.atoms()[i_at];
              if (atom.uij_is_defined()) n_anisou++;
              model_atom_labels.push_back(atom.pdb_label_columns_segid_small_str());
              i_model_atom++;
              element_charge_types[atom.pdb_element_charge_columns()]++;
            }
          }
          {
            typedef std::set<std::string>::const_iterator it;
            it i_end = rg_resnames.end();
            for(it i=rg_resnames.begin();i!=i_end;i++) {
              resnames[std::string((*i).c_str())]++;
            }
          }
          {
            typedef std::map<std::string, std::vector<std::string> >::const_iterator it;
            it i_end = altloc_resnames.end();
            for(it i=altloc_resnames.begin();i!=i_end;i++) {
              if (i->second.size() != 1U) {
                residue_groups_with_multiple_resnames_using_same_altloc
                  .push_back(rg);
              }
            }
            // This should do the same as above, C++11 way.
            // for (auto i: altloc_resnames) {
            //   if (i->second.size() != 1U) {
            //     residue_groups_with_multiple_resnames_using_same_altloc
            //       .push_back(rg);
            //   }
            // }
          }
          if (have_blank_altloc) {
            n_alt_conf_improper++;
            if (!chain_alt_conf_improper) {
              chain_alt_conf_improper = boost::optional<residue_group>(rg);
            }
          }
          else if (have_main_conf) {
            if (rg_altlocs.size() == 0) {
              n_alt_conf_none++;
            }
            else {
              n_alt_conf_proper++;
              if (!chain_alt_conf_proper) {
                chain_alt_conf_proper = boost::optional<residue_group>(rg);
              }
            }
          }
          else if (rg_altlocs.size() != 0) {
            n_alt_conf_pure++;
          }
          chain_altlocs.insert(rg_altlocs.begin(), rg_altlocs.end());
          if (rg_resnames.size() == 1) {
            n_residues++;
          }
          else if (rg_resnames.size() != 0) {
            n_residue_groups++;
          }
          if (prev_rg && prev_rg->resid() == rg.resid()) {
            consecutive_residue_groups_with_same_resid.push_back(
              af::tiny<residue_group, 2>(*prev_rg, rg));
          }
          prev_rg = boost::optional<residue_group>(rg);
        }
        {
          typedef std::set<std::string>::const_iterator it;
          it i_end = chain_altlocs.end();
          for(it i=chain_altlocs.begin();i!=i_end;i++) {
            alt_conf_ids[*i]++;
            n_alt_conf++;
          }
        }
        if (   chain_alt_conf_proper
            && chain_alt_conf_improper) {
          n_chains_with_mix_of_proper_and_improper_alt_conf++;
          if (   !alt_conf_proper
              || !alt_conf_improper) {
            alt_conf_proper = chain_alt_conf_proper;
            alt_conf_improper = chain_alt_conf_improper;
          }
        }
        else {
          if (!alt_conf_proper) alt_conf_proper = chain_alt_conf_proper;
          if (!alt_conf_improper) alt_conf_improper = chain_alt_conf_improper;
        }
      }
      {
        typedef std::map<std::string, unsigned>::const_iterator it;
        it i_end = model_chain_ids.end();
        for(it i=model_chain_ids.begin();i!=i_end;i++) {
          unsigned count = i->second;
          if (count != 1) n_duplicate_chain_ids += count;
        }
      }
      IOTBX_ASSERT(i_model_atom == model_atoms_size);
      n_atoms += model_atoms_size;
      n_duplicate_atom_labels += detail::find_duplicate_atom_labels(
        duplicate_atom_labels,
        model,
        model_atoms_size,
        model_atom_labels);
    }
    {
      typedef std::map<std::string, unsigned>::const_iterator it;
      it i_end = model_ids.end();
      for(it i=model_ids.begin();i!=i_end;i++) {
        unsigned count = i->second;
        if (count != 1) n_duplicate_model_ids += count;
      }
    }
    {
      typedef std::map<std::string, unsigned>::const_iterator it;
      it i_end = resnames.end();
      for(it i=resnames.begin();i!=i_end;i++) {
        resname_classes[common_residue_names::get_class(i->first)] +=i->second;
      }
    }
  }

}}} // namespace iotbx::pdb::hierarchy
