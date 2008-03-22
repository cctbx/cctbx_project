#include <iotbx/pdb/hierarchy_v2.h>

namespace iotbx { namespace pdb { namespace hierarchy_v2 {

namespace detail {

  void
  append_range(
    std::vector<unsigned>& v,
    unsigned i_seq,
    unsigned i_seq_end)
  {
    for(;i_seq!=i_seq_end;i_seq++) v.push_back(i_seq);
  }

} // namespace detail

  atom_selection_cache::atom_selection_cache(
    hierarchy_v2::root const& root)
  {
    unsigned i_seq = 0;
    std::vector<model> const& models = root.models();
    unsigned n_mds = root.models_size();
    for(unsigned i_md=0;i_md<n_mds;i_md++) {
      unsigned model_i_seq_start = i_seq;
      model const& mo = models[i_md];
      unsigned n_chs = mo.chains_size();
      std::vector<chain> const& chains = mo.chains();
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) {
        unsigned chain_i_seq_start = i_seq;
        chain const& ch = chains[i_ch];
        unsigned n_rgs = ch.residue_groups_size();
        std::vector<residue_group> const& rgs = ch.residue_groups();
        for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) {
          unsigned rg_i_seq_start = i_seq;
          residue_group const& rg = rgs[i_rg];
          unsigned n_ags = rg.atom_groups_size();
          std::vector<atom_group> const& ags = rg.atom_groups();
          for(unsigned i_ag=0;i_ag<n_ags;i_ag++) {
            unsigned ag_i_seq_start = i_seq;
            atom_group const& ag = ags[i_ag];
            unsigned n_ats = ag.atoms_size();
            std::vector<atom> const& ats = ag.atoms();
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              name[ats[i_at].data->name.elems].push_back(i_seq++);
            }
            i_seq = ag_i_seq_start;
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              segid[ats[i_at].data->segid.elems].push_back(i_seq++);
            }
            i_seq = ag_i_seq_start;
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              element[ats[i_at].data->element.elems].push_back(i_seq++);
            }
            i_seq = ag_i_seq_start;
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              charge[ats[i_at].data->charge.elems].push_back(i_seq++);
            }
            i_seq = ag_i_seq_start;
            for(unsigned i_at=0;i_at<n_ats;i_at++,i_seq++) {
              if (ats[i_at].uij_is_defined()) anisou.push_back(i_seq);
            }
            detail::append_range(
              altloc[ag.data->altloc.elems[0] == '\0'
                ? " " : ag.data->altloc.elems], ag_i_seq_start, i_seq);
            detail::append_range(
              resname[ag.data->resname.elems], ag_i_seq_start, i_seq);
          }
          detail::append_range(
            resseq[rg.data->resseq.elems], rg_i_seq_start, i_seq);
          detail::append_range(
            icode[rg.data->icode.elems], rg_i_seq_start, i_seq);
          detail::append_range(
            resid[rg.resid()], rg_i_seq_start, i_seq);
        }
        detail::append_range(chain_id[ch.data->id], chain_i_seq_start, i_seq);
      }
      detail::append_range(model_id[mo.data->id], model_i_seq_start, i_seq);
    }
    n_seq = i_seq;
  }

}}} // namespace iotbx::pdb::hierarchy_v2
