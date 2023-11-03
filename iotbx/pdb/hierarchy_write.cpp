#include <iotbx/pdb/hierarchy_atoms.h>
#include <iotbx/pdb/write_utils.h>

namespace iotbx { namespace pdb { namespace hierarchy {

  void
  residue_groups_as_pdb_string(
    stream_write& write,
    atom_label_columns_formatter& label_formatter,
    std::vector<residue_group> const& residue_groups,
    int interleaved_conf,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij,
    bool output_break_records)
  {
    char buf[81U + 81U + 81U + 81U];
    unsigned n_rg = static_cast<unsigned>(residue_groups.size());
    for(unsigned i_rg=0;i_rg<n_rg;i_rg++) {
      residue_group const& rg = residue_groups[i_rg];
      if (i_rg != 0 && !rg.data->link_to_previous && output_break_records) {
        write("BREAK\n", 6U);
      }
      label_formatter.resseq = rg.data->resseq.elems;
      label_formatter.icode = rg.data->icode.elems;
      if (interleaved_conf <= 0) {
        unsigned n_ag = rg.atom_groups_size();
        for(unsigned i_ag=0;i_ag<n_ag;i_ag++) {
          atom_group const& ag = rg.atom_groups()[i_ag];
          label_formatter.altloc = ag.data->altloc.c_str();
          label_formatter.resname = ag.data->resname.c_str();
          typedef std::vector<atom> va;
          va const& atoms = ag.atoms();
          va::const_iterator atoms_end = atoms.end();
          for(va::const_iterator
                atom=atoms.begin();atom!=atoms_end;atom++) {
#define IOTBX_LOC \
            unsigned str_len = atom->format_atom_record_group( \
              buf, &label_formatter, \
              atom_hetatm, sigatm, anisou, siguij); \
            if (str_len != 0) { \
              buf[str_len] = '\n'; \
              write(buf, str_len+1U); \
            }
            IOTBX_LOC
          }
        }
      }
      else { // interleaved_conf > 0
        af::shared<atom> atoms_ilc = rg.atoms_interleaved_conf(
          /* group_residue_names */ (interleaved_conf < 2));
        af::const_ref<atom> ats = atoms_ilc.const_ref();
        unsigned n_at = static_cast<unsigned>(ats.size());
        for(unsigned i_at=0;i_at<n_at;i_at++) {
          hierarchy::atom const* atom = &ats[i_at];
          shared_ptr<atom_group_data> ag_data = atom->parent_ptr();
          label_formatter.altloc = ag_data->altloc.c_str();
          label_formatter.resname = ag_data->resname.c_str();
          IOTBX_LOC
        }
      }
#undef IOTBX_LOC
    }
  }

  // Why we have separate write functions? see
  // iotbx/pdb/construct_hierarchy.cpp:input_atoms_with_labels_generator::run
  void
  models_as_pdb_string(
    stream_write& write,
    std::vector<model> const& models,
    bool append_end,
    int interleaved_conf,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij,
    bool output_break_records)
  {
    atom_label_columns_formatter label_formatter;
    unsigned n_mds = static_cast<unsigned>(models.size());
    for(unsigned i_md=0;i_md<n_mds;i_md++) {
      if (n_mds != 1U) {
        write_utils::model_record(write, models[i_md].data->id);
      }
      unsigned n_chs = models[i_md].chains_size();
      std::vector<chain> const& chains = models[i_md].chains();
      for(unsigned i_ch=0;i_ch<n_chs;i_ch++) {
        chain const& ch = chains[i_ch];
        label_formatter.chain_id = ch.data->id.c_str();
        residue_groups_as_pdb_string(
          write,
          label_formatter,
          ch.residue_groups(),
          interleaved_conf,
          atom_hetatm, sigatm, anisou, siguij, output_break_records);
        if (ch.is_polymer_chain()) write("TER\n", 4U);
      }
      if (n_mds != 1U) {
        write("ENDMDL\n", 7U);
      }
    }
    if (append_end) {
      write("END\n", 4U);
    }
  }

  void
  root::write_pdb_file(
    const char* file_name,
    bool open_append,
    bool append_end,
    int interleaved_conf,
    boost::optional<int> const& atoms_reset_serial_first_value,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij,
    bool output_break_records) const
  {
    if (atoms_reset_serial_first_value) {
      atoms_reset_serial(interleaved_conf, *atoms_reset_serial_first_value);
    }
    write_utils::fstream_open_close foc(file_name, open_append);
    write_utils::fstream_write write(&foc.out);
    models_as_pdb_string(
      write,
      models(),
      append_end,
      interleaved_conf,
      atom_hetatm,
      sigatm,
      anisou,
      siguij,
      output_break_records);
  }

}}} // namespace iotbx::pdb::hierarchy
