#include <iotbx/pdb/hierarchy_v2.h>
#include <iotbx/pdb/common_residue_names.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 {

namespace detail {

  bool
  cmp_first_elements_of_vectors(
    const std::vector<unsigned>* a,
    const std::vector<unsigned>* b)
  {
    return ((*a)[0] < (*b)[0]);
  }

} // namespace detail

  overall_counts::overall_counts(
    hierarchy_v2::root const& root_)
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
    n_alt_conf_none(0),
    n_alt_conf_pure(0),
    n_alt_conf_proper(0),
    n_alt_conf_improper(0),
    n_chains_with_mix_of_proper_and_improper_alt_conf(0)
  {
    af::shared<atom> atoms_owner = root.atoms();
    n_atoms = static_cast<unsigned>(atoms_owner.size());
    const atom* atoms = atoms_owner.begin();
    atoms_reset_tmp(atoms_owner.const_ref());
    for(unsigned i_md=0;i_md<n_models;i_md++) {
      hierarchy_v2::model const& model = root.models()[i_md];
      if (model.chains_size() == 0) n_empty_models++;
      model_ids[model.data->id]++;
      std::map<std::string, unsigned> model_chain_ids;
      std::map<std::string, std::vector<unsigned> > model_atom_labels_i_seqs;
      unsigned n_ch = model.chains_size();
      n_chains += n_ch;
      for(unsigned i_ch=0;i_ch<n_ch;i_ch++) {
        hierarchy_v2::chain const& chain = model.chains()[i_ch];
        if (chain.residue_groups_size() == 0) n_empty_chains += 1;
        model_chain_ids[chain.data->id]++;
        chain_ids[chain.data->id]++;
        std::set<char> chain_altlocs;
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
          std::set<char> rg_altlocs;
          std::set<str3> rg_resnames;
          unsigned n_ag = rg.atom_groups_size();
          for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
            atom_group const& ag = rg.atom_groups()[i_ag];
            if (ag.atoms_size() == 0) n_empty_atom_groups++;
            char altloc = ag.data->altloc.elems[0];
            if (altloc == '\0') {
              have_main_conf = true;
            }
            else {
              if (altloc == blank_altloc_char) {
                have_blank_altloc = true;
              }
              rg_altlocs.insert(altloc);
            }
            rg_resnames.insert(ag.data->resname);
            unsigned n_atoms = ag.atoms_size();
            for(unsigned i_at=0;i_at<n_atoms;i_at++) {
              hierarchy_v2::atom const& atom = ag.atoms()[i_at];
              model_atom_labels_i_seqs[atom.pdb_label_columns()]
                .push_back(atom.data->tmp);
              element_charge_types[atom.pdb_element_charge_columns()]++;
            }
          }
          {
            typedef std::set<str3>::const_iterator it;
            it i_end = rg_resnames.end();
            for(it i=rg_resnames.begin();i!=i_end;i++) {
              resnames[std::string((*i).elems)]++;
            }
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
          typedef std::set<char>::const_iterator it;
          it i_end = chain_altlocs.end();
          for(it i=chain_altlocs.begin();i!=i_end;i++) {
            alt_conf_ids[std::string(&*i, 1U)]++;
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
      {
        std::vector<std::vector<unsigned>*> model_duplicate_atom_labels;
        typedef std::map<std::string, std::vector<unsigned> >::iterator it;
        it i_end = model_atom_labels_i_seqs.end();
        for(it i=model_atom_labels_i_seqs.begin();i!=i_end;i++) {
          if (i->second.size() == 1) continue;
          model_duplicate_atom_labels.push_back(&(i->second));
        }
        std::size_t n_mdal = model_duplicate_atom_labels.size();
        if (n_mdal != 0) {
          std::sort(
            model_duplicate_atom_labels.begin(),
            model_duplicate_atom_labels.end(),
            detail::cmp_first_elements_of_vectors);
          duplicate_atom_labels.reserve(duplicate_atom_labels.size() + n_mdal);
          for(std::size_t i_mdal=0;i_mdal<n_mdal;i_mdal++) {
            std::vector<unsigned>&
              i_seqs = *model_duplicate_atom_labels[i_mdal];
            unsigned n_is = static_cast<unsigned>(i_seqs.size());
            n_duplicate_atom_labels += n_is;
            af::shared<atom> group((af::reserve(n_is)));
            for(unsigned i_is=0;i_is<n_is;i_is++) {
              group.push_back(atoms[i_seqs[i_is]]);
            }
            duplicate_atom_labels.push_back(group);
          }
        }
      }
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

}}} // namespace iotbx::pdb::hierarchy_v2
