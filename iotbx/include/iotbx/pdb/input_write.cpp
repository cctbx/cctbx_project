#include <iotbx/pdb/input.h>
#include <iotbx/pdb/write_utils.h>

namespace iotbx { namespace pdb {

  void
  input_as_pdb_string(
    stream_write& write,
    input const& self,
    bool append_end,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij)
  {
    af::const_ref<std::string>
      model_ids = self.model_ids().const_ref();
    af::const_ref<std::vector<unsigned> >
      chain_indices = self.chain_indices().const_ref();
    SCITBX_ASSERT(chain_indices.size() == model_ids.size());
    const std::size_t* break_index = self.break_indices().begin();
    const std::size_t* break_indices_end = self.break_indices().end();
    unsigned next_break_index = static_cast<unsigned>(
      break_index == break_indices_end ?
        self.atoms().size() : *break_index++);
    const input_atom_labels* iall = self.input_atom_labels_list().begin();
    const hierarchy::atom* atoms = self.atoms().begin();
    unsigned next_chain_range_begin = 0;
    for(unsigned i_model=0;i_model<model_ids.size();i_model++) {
      std::string const& model_id = model_ids[i_model];
      if (model_id.size() != 0) {
        write_utils::model_record(write, model_id);
      }
      range_loop<unsigned> ch_r(
        chain_indices[i_model], next_chain_range_begin);
      for(unsigned i_chain=0;ch_r.next();i_chain++) {
        bool is_first_in_chain = true;
        for (unsigned i_atom=ch_r.begin; i_atom!=ch_r.end; i_atom++) {
          input_atom_labels const& ial = iall[i_atom];
          bool is_first_after_break = (i_atom == next_break_index);
          if (is_first_after_break) {
            next_break_index = static_cast<unsigned>(
              break_index == break_indices_end ?
                self.atoms().size() : *break_index++);
            write("BREAK\n", 6U);
          }
          hierarchy::atom_with_labels awl(
            atoms[i_atom],
            model_id.c_str(),
            ial.chain_small().elems,
            ial.resseq_small().elems,
            ial.icode_small().elems,
            ial.altloc_small().elems,
            ial.resname_small().elems,
            is_first_in_chain,
            is_first_after_break);
          is_first_in_chain = false;
          std::string s = awl.format_atom_record_group(
            atom_hetatm, sigatm, anisou, siguij);
          write(s.c_str(), s.size());
          write("\n", 1U);
        }
        write("TER\n", 4U);
      }
      if (model_id.size() != 0) {
        write("ENDMDL\n", 7U);
      }
      next_chain_range_begin = ch_r.end;
    }
    SCITBX_ASSERT(break_index == break_indices_end);
    if (append_end) {
      write("END\n", 4U);
    }
  }

  void
  input::write_pdb_file(
    const char* file_name,
    bool open_append,
    bool append_end,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij) const
  {
    write_utils::fstream_open_close foc(file_name, open_append);
    write_utils::fstream_write write(&foc.out);
    input_as_pdb_string(
      write, *this, append_end, atom_hetatm, sigatm, anisou, siguij);
  }

}} // namespace iotbx::pdb
