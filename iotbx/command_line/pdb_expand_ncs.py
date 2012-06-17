# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb.expand_ncs

import libtbx.phil
from libtbx.utils import Sorry, Usage
import os
import sys

master_phil = libtbx.phil.parse("""
expand_ncs {
  file_name = None
    .type = path
    .short_caption = PDB file
    .style = file_type:pdb bold
  output_file = None
    .type = path
    .short_caption = Output file
    .style = file_type:pdb new_file
  clash_limit = 20
    .type = int
}
""")

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  if (params is None) :
    if (len(args) == 0) :
      raise Usage("iotbx.pdb.expand_ncs model.pdb")
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="expand_ncs.file_name")
    params = cmdline.work.extract()
  from iotbx import file_reader
  import iotbx.pdb.hierarchy
  pdb_file = file_reader.any_file(params.expand_ncs.file_name,
    force_type="pdb")
  pdb_in = pdb_file.file_object
  matrices = pdb_in.process_mtrix_records()
  pdb_hierarchy = pdb_in.construct_hierarchy()
  hierarchy_new = iotbx.pdb.hierarchy.expand_ncs(
    pdb_hierarchy=pdb_hierarchy,
    matrices=matrices,
    log=out)
  if (params.expand_ncs.output_file is None) :
    params.expand_ncs.output_file = os.path.basename(
      os.path.splitext(params.expand_ncs.file_name)[0]) + "_ncs.pdb"
  f = open(params.expand_ncs.output_file, "w")
  f.write(hierarchy_new.as_pdb_string(pdb_in.crystal_symmetry()))
  f.close()
  print >> out, "Wrote %s" % params.expand_ncs.output_file
  if (params.expand_ncs.clash_limit is not None) :
    n_bad_pairs = iotbx.pdb.quick_clash_check(
      file_name=params.expand_ncs.output_file)
    if (n_bad_pairs > params.expand_ncs.clash_limit) :
      raise Sorry(("%d bad clashes found in output model (%s).  Set "+
        "clash_limit=None to disable the clash check.") %
        (n_bad_pairs, params.expand_ncs.output_file))
  return os.path.abspath(params.expand_ncs.output_file)

def validate_params (params) :
  if (params.expand_ncs.file_name is None) :
    raise Sorry("No PDB file specified.")
  elif (not os.path.isfile(params.expand_ncs.file_name)) :
    raise Sorry("%s is not a valid file." % params.expand_ncs.file_name)
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
