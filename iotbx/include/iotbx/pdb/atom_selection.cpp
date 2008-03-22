#include <iotbx/pdb/hierarchy_v2.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 {

  atom_selection_cache::atom_selection_cache(
    hierarchy_v2::root const& root)
  {
    unsigned i_seq = 0;
    std::vector<model> const& models = root.models();
    unsigned n_mds = root.models_size();
    for(unsigned i_md=0;i_md<n_mds;i_md++) {
      model const& mo = models[i_md];
      unsigned n_chs = mo.chains_size();
      std::vector<chain> const& chains = mo.chains();
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) {
        chain const& ch = chains[i_ch];
        unsigned n_rgs = ch.residue_groups_size();
        std::vector<residue_group> const& rgs = ch.residue_groups();
        for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) {
          residue_group const& rg = rgs[i_rg];
          std::string rg_resseq = rg.data->resseq.elems;
          std::string rg_icode = rg.data->icode.elems;
          std::string rg_resid = rg.resid();
          unsigned n_ags = rg.atom_groups_size();
          std::vector<atom_group> const& ags = rg.atom_groups();
          for(unsigned i_ag=0;i_ag<n_ags;i_ag++) {
            atom_group const& ag = ags[i_ag];
            std::string ag_altloc = (ag.data->altloc.elems[0] == '\0'
              ? " " : ag.data->altloc.elems);
            std::string ag_resname = ag.data->resname.elems;
            unsigned n_ats = ag.atoms_size();
            std::vector<atom> const& ats = ag.atoms();
            for(unsigned i_at=0;i_at<n_ats;i_at++,i_seq++) {
              hierarchy_v2::atom const& atom = ats[i_at];
              atom_data const& ad = *atom.data;
              name[ad.name.elems].push_back(i_seq);
              altloc[ag_altloc].push_back(i_seq);
              resname[ag_resname].push_back(i_seq);
              chain_id[ch.data->id].push_back(i_seq);
              resseq[rg_resseq].push_back(i_seq);
              icode[rg_icode].push_back(i_seq);
              resid[rg_resid].push_back(i_seq);
              segid[ad.segid.elems].push_back(i_seq);
              model_id[mo.data->id].push_back(i_seq);
              element[ad.element.elems].push_back(i_seq);
              charge[ad.charge.elems].push_back(i_seq);
              if (atom.uij_is_defined()) anisou.push_back(i_seq);
            }
          }
        }
      }
    }
    n_seq = i_seq;
  }

}}} // namespace iotbx::pdb::hierarchy_v2
