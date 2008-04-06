#include <iotbx/pdb/hierarchy.h>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb { namespace hierarchy {

  root
  root::select(
    af::const_ref<bool> const& atom_selection) const
  {
    root result;
    unsigned n_sel = static_cast<unsigned>(atom_selection.size());
    unsigned i_sel = 0;
    unsigned n_mds = models_size();
    for(unsigned i_md=0;i_md<n_mds;i_md++) {
      model const& md = data->models[i_md];
      model r_md(md.data->id);
      unsigned n_chs = md.chains_size();
      std::vector<chain> const& chs = md.chains();
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) {
        chain const& ch = chs[i_ch];
        chain r_ch(ch.data->id);
        unsigned n_rgs = ch.residue_groups_size();
        std::vector<residue_group> const& rgs = ch.residue_groups();
        for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) {
          residue_group const& rg = rgs[i_rg];
          residue_group r_rg(rg.data->resseq.elems, rg.data->icode.elems);
          unsigned n_ags = rg.atom_groups_size();
          std::vector<atom_group> const& ags = rg.atom_groups();
          for(unsigned i_ag=0;i_ag<n_ags;i_ag++) {
            atom_group const& ag = ags[i_ag];
            unsigned n_ats = ag.atoms_size();
            if (i_sel + n_ats > n_sel) {
              throw std::invalid_argument("atom_selection array too short.");
            }
            atom_group r_ag(ag.data->altloc.elems, ag.data->resname.elems);
            std::vector<atom> const& atoms = ag.atoms();
            boost::scoped_array<atom> r_atoms(new atom[n_ats]);
            unsigned r_n_ats = 0;
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              if (atom_selection[i_sel++]) r_atoms[r_n_ats++] = atoms[i_at];
            }
            if (r_n_ats != 0) {
              r_ag.pre_allocate_atoms(r_n_ats);
              for(unsigned i_at=0;i_at<r_n_ats;i_at++) {
                atom r_at = r_atoms[i_at].detached_copy();
                r_ag.append_atom(r_at);
              }
              r_rg.append_atom_group(r_ag);
            }
          }
          if (r_rg.atom_groups_size() != 0) {
            r_ch.append_residue_group(r_rg);
          }
        }
        if (r_ch.residue_groups_size() != 0) {
          r_md.append_chain(r_ch);
        }
      }
      if (r_md.chains_size() != 0) {
        result.append_model(r_md);
      }
    }
    if (i_sel != n_sel) {
      throw std::invalid_argument("atom_selection array too large.");
    }
    return result;
  }

}}} // namespace iotbx::pdb::hierarchy
