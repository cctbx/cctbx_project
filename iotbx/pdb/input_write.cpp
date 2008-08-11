#include <iotbx/pdb/input.h>
#include <iotbx/pdb/write_utils.h>

namespace iotbx { namespace pdb {

namespace {

  struct input_as_pdb_string_generator : input_atoms_with_labels_generator
  {
    stream_write& write;
    bool append_end;
    bool atom_hetatm;
    bool sigatm;
    bool anisou;
    bool siguij;

    public:
      input_as_pdb_string_generator(
        stream_write& write_,
        bool append_end_,
        bool atom_hetatm_,
        bool sigatm_,
        bool anisou_,
        bool siguij_)
      :
        write(write_),
        append_end(append_end_),
        atom_hetatm(atom_hetatm_),
        sigatm(sigatm_),
        anisou(anisou_),
        siguij(siguij_)
      {}

      virtual bool
      process_model(str8 const& model_id)
      {
        if (model_id.size() != 0) {
          write_utils::model_record(write, model_id);
        }
        return true;
      }

      virtual bool
      process_endmdl(str8 const& model_id)
      {
        if (model_id.size() != 0) {
          write("ENDMDL\n", 7U);
        }
        return true;
      }

      virtual bool
      process_atom(hierarchy::atom_with_labels const& awl)
      {
        std::string s = awl.format_atom_record_group(
          atom_hetatm, sigatm, anisou, siguij);
        write(s.c_str(), s.size());
        write("\n", 1U);
        return true;
      }

      virtual bool
      process_break()
      {
        write("BREAK\n", 6U);
        return true;
      }

      virtual bool
      process_ter()
      {
        write("TER\n", 4U);
        return true;
      }

      virtual void
      process_end()
      {
        if (append_end) {
          write("END\n", 4U);
        }
      }
  };

} // namespace <anonymous>

  void
  input_as_pdb_string(
    input const& self,
    stream_write& write,
    bool append_end,
    bool atom_hetatm,
    bool sigatm,
    bool anisou,
    bool siguij)
  {
    input_as_pdb_string_generator(
      write, append_end, atom_hetatm, sigatm, anisou, siguij).run(self);
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
    input_as_pdb_string_generator(
      write, append_end, atom_hetatm, sigatm, anisou, siguij).run(*this);
  }

}} // namespace iotbx::pdb
