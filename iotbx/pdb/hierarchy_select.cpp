#include <iotbx/pdb/hierarchy.h>
#include <boost/scoped_array.hpp>

namespace iotbx { namespace pdb { namespace hierarchy {

  root
  root::select(
    af::const_ref<bool> const& atom_selection,
    bool copy_atoms) const
  {
    root result;
    unsigned n_sel = static_cast<unsigned>(atom_selection.size());
    unsigned i_sel = 0;
#define IOTBX_LOC_OPEN \
    unsigned n_mds = models_size(); \
    for(unsigned i_md=0;i_md<n_mds;i_md++) { \
      model const& md = data->models[i_md]; \
      model r_md(md.data->id.c_str()); \
      unsigned n_chs = md.chains_size(); \
      std::vector<chain> const& chs = md.chains(); \
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) { \
        chain const& ch = chs[i_ch]; \
        chain r_ch(ch.data->id.c_str()); \
        unsigned n_rgs = ch.residue_groups_size(); \
        std::vector<residue_group> const& rgs = ch.residue_groups(); \
        for(unsigned i_rg=0;i_rg<n_rgs;i_rg++) { \
          residue_group const& rg = rgs[i_rg]; \
          residue_group r_rg(rg.data->resseq.elems, rg.data->icode.elems); \
          unsigned n_ags = rg.atom_groups_size(); \
          std::vector<atom_group> const& ags = rg.atom_groups(); \
          for(unsigned i_ag=0;i_ag<n_ags;i_ag++) { \
            atom_group const& ag = ags[i_ag]; \
            unsigned n_ats = ag.atoms_size(); \
            atom_group r_ag(ag.data->altloc, ag.data->resname); \
            std::vector<atom> const& atoms = ag.atoms(); \
            boost::scoped_array<atom> r_atoms(new atom[n_ats]); \
            unsigned r_n_ats = 0;
            IOTBX_LOC_OPEN
            if (i_sel + n_ats > n_sel) {
              throw std::invalid_argument("atom_selection array too short.");
            }
            for(unsigned i_at=0;i_at<n_ats;i_at++) {
              if (atom_selection[i_sel++]) r_atoms[r_n_ats++] = atoms[i_at];
            }
#define IOTBX_LOC_CLOSE(BREAK) \
            if (r_n_ats != 0) { \
              r_ag.pre_allocate_atoms(r_n_ats); \
              for(unsigned i_at=0;i_at<r_n_ats;i_at++) { \
                if (copy_atoms) { \
                  atom r_at = r_atoms[i_at].detached_copy(); \
                  r_ag.append_atom(r_at); \
                } \
                else { \
                  r_ag.append_atom_with_other_parent(r_atoms[i_at]); \
                } \
              } \
              r_rg.append_atom_group(r_ag); \
              BREAK \
            } \
          } \
          if (r_rg.atom_groups_size() != 0) { \
            r_ch.append_residue_group(r_rg); \
            BREAK \
          } \
        } \
        if (r_ch.residue_groups_size() != 0) { \
          r_md.append_chain(r_ch); \
          BREAK \
        } \
      } \
      if (r_md.chains_size() != 0) { \
        result.append_model(r_md); \
        BREAK \
      } \
    }
#define IOTBX_LOC_BREAK
    IOTBX_LOC_CLOSE(IOTBX_LOC_BREAK)
#undef IOTBX_LOC_BREAK
    if (i_sel != n_sel) {
      throw std::invalid_argument("atom_selection array too large.");
    }
    return result;
  }

  root
  root::select(
    af::const_ref<std::size_t> const& atom_selection,
    bool copy_atoms) const
  {
    root result;
    const std::size_t* a_s_end = atom_selection.end();
    const std::size_t* a_s = atom_selection.begin();
    if (a_s == a_s_end) return result;
    std::size_t i_seq = 0;
    IOTBX_LOC_OPEN
      for(unsigned i_at=0;i_at<n_ats;i_at++,i_seq++) {
        if (*a_s == i_seq) {
          r_atoms[r_n_ats++] = atoms[i_at];
          if (++a_s == a_s_end) break;
          if (*a_s <= i_seq) {
            throw std::invalid_argument(
              "atom_selection indices not in strictly ascending order.");
          }
        }
      }
#define IOTBX_LOC_BREAK if (a_s == a_s_end) break;
    IOTBX_LOC_CLOSE(IOTBX_LOC_BREAK)
#undef IOTBX_LOC_BREAK
    if (a_s != a_s_end) {
      throw std::invalid_argument(
        "atom_selection indices greater than or equal to number of atoms.");
    }
    return result;
  }

#undef IOTBX_LOC_OPEN
#undef IOTBX_LOC_CLOSE

}}} // namespace iotbx::pdb::hierarchy
