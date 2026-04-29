"""Combine models into a single model"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.cif.model
import iotbx.phil
import libtbx
from libtbx.utils import Usage, format_cpu_times
import sys, os

master_phil = iotbx.phil.parse("""
join_fragment_files {
  reset_atom_serial = True
    .type = bool
  model_file = None
    .type = path
    .multiple = True
  format = mmcif pdb
    .type = choice
}
""")

def run(args,  command_name="iotbx.pdb.join_fragment_files"):
  from iotbx import file_reader
  def usage():
    raise Usage("""\
%s file1.pdb file2.pdb [...]

or define the environment variable
  PDB_MIRROR_PDB
to join all fragment files in the PDB.""" % command_name)

  if (len(args) == 0 or args == ["--exercise"]):
    pdb_mirror_pdb = os.environ.get("PDB_MIRROR_PDB")
    if (pdb_mirror_pdb is None):
      if (len(args) == 0): usage()
    else:
      for line in iotbx.pdb.pdb_codes_fragment_files.splitlines():
        print("PDB code group:", line)
        codes = line.split()
        joined = iotbx.pdb.join_fragment_files(
          file_names = [
            os.path.join(pdb_mirror_pdb, code[1:3], "pdb%s.ent.gz" % code)
              for code in codes]).joined
        file_name_out = "%s_%s.pdb" % (codes[0], codes[-1])
        print("  writing:", file_name_out)
        out = open(file_name_out, "w")
        print("\n".join(joined.info), file=out)
        out.write(joined.as_pdb_string(append_end=True))
        if (len(args) != 0): break
    print(format_cpu_times())
  else:
    sources = []
    file_names = []
    interpreter = master_phil.command_line_argument_interpreter()
    input_file_type = None
    for arg in args :
      if os.path.isfile(arg):
        input_file = file_reader.any_file(arg)
        if (input_file.file_type == "pdb"):
          file_names.append(input_file)
          sources.append(interpreter.process(arg="model_file=\"%s\"" % arg))
        elif (input_file.file_type == "phil"):
          sources.append(input_file.file_object)
      else :
        arg_phil = interpreter.process(arg=arg)
        sources.append(arg_phil)
    work_phil = master_phil.fetch(sources=sources)
    work_params = work_phil.extract()
    file_names = work_params.join_fragment_files.model_file
    if (len(file_names) < 2): usage()
    result = iotbx.pdb.join_fragment_files(file_names=file_names)
    joined = result.joined
    if work_params.join_fragment_files.reset_atom_serial:
      joined.atoms_reset_serial()
    if work_params.join_fragment_files.format in (None, libtbx.Auto, "pdb"):
      print("\n".join(joined.info))
      sys.stdout.write(joined.as_pdb_string(append_end=True))
    elif work_params.join_fragment_files.format == "mmcif":
      cif_object = iotbx.cif.model.cif()
      cif_object["combined"] = joined.as_cif_block(
        crystal_symmetry=result.crystal_symmetry)
      cif_object.show(out=sys.stdout)

if (__name__ == "__main__"):
  run(sys.argv[1:])
