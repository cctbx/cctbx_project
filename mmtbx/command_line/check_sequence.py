
import mmtbx.validation.sequence
from libtbx.utils import Sorry, Usage
import libtbx.phil
import os
import sys

master_phil = libtbx.phil.parse("""
input
  .style = auto_align
{
  pdb_file = None
    .type = path
    .style = file_type:pdb input_file
  seq_file = None
    .type = path
    .style = file_type:seq input_file
}
include scope mmtbx.validation.sequence.master_phil
""", process_includes=True)

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  if (len(args) == 0) and (params is None) :
    raise Usage("""\
mmtbx.check_sequence model.pdb sequence.fa

Verify the sequence of each chain in a PDB file to detect residue mismatches
and other inconsistencies (similar to validation upon PDB deposition).""")
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="input.pdb_file",
    seq_file_def="input.seq_file")
  params = cmdline.work.extract()
  try :
    validate_params(params)
  except Sorry, e :
    print e
    raise Usage("mmtbx.check_sequence model.pdb sequence.fa")
  import mmtbx.validation.sequence
  from iotbx.file_reader import any_file
  pdb_in = any_file(params.input.pdb_file, force_type="pdb")
  pdb_in.check_file_type("pdb")
  seq_in = any_file(params.input.seq_file, force_type="seq")
  seq_in.check_file_type("seq")
  pdb_hierarchy = pdb_in.file_object.construct_hierarchy()
  sequences = seq_in.file_object
  if (len(sequences) == 0) :
    raise Sorry("There don't appear to be any valid sequences in %s!" %
      params.input.seq_file)
  v = mmtbx.validation.sequence.validation(
    pdb_hierarchy=pdb_hierarchy,
    sequences=sequences,
    params=params,
    log=out)
  v.show(out=out)

def validate_params (params) :
  if (params.input.pdb_file is None) :
    raise Sorry("No PDB file specified.")
  elif (not os.path.isfile(params.input.pdb_file)) :
    raise Sorry("'%s' is not a file." % params.input.pdb_file)
  elif (params.input.seq_file is None) :
    raise Sorry("No sequence file specified.")
  elif (not os.path.isfile(params.input.seq_file)) :
    raise Sorry("'%s' is not a file." % params.input.seq_file)
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
