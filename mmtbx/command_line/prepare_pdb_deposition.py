from __future__ import division
from __future__ import print_function
import os
import sys

import iotbx.phil
from iotbx.file_reader import any_file
from iotbx.pdb import mmcif
from mmtbx.command_line import model_vs_sequence
from iotbx.cif import category_sort_function


master_phil = iotbx.phil.parse("""
include scope mmtbx.command_line.cc_star.master_phil
""", process_includes=True)

master_phil = iotbx.phil.parse("""
input
  .style = auto_align
{
  pdb_file = None
    .type = path
    .style = file_type:pdb input_file
  seq_file = None
    .type = path
    .style = file_type:seq input_file
  include scope mmtbx.command_line.cc_star.master_phil
}
output
  .style = auto_align
{
  cif_file = None
    .type = path
}
high_resolution = None
  .type = float
  .input_size = 64
low_resolution = None
  .type = float
  .input_size = 64
include scope mmtbx.validation.sequence.master_phil
""", process_includes=True)


def run(args, out=None):
  import iotbx.phil
  if (out is None) :
    out = sys.stdout
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="input.pdb_file",
    seq_file_def="input.seq_file"
  )
  params = cmdline.work.extract()
  cmdline.work.show()
  if params.output.cif_file is None:
    params.output.cif_file = os.path.splitext(params.input.pdb_file)[0] + ".deposit.cif"
  model_vs_sequence.validate_params(params)
  pdb_input = mmcif.cif_input(file_name=params.input.pdb_file)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  cif_model = pdb_input.cif_model
  cif_block = pdb_input.cif_block
  seq_in = any_file(params.input.seq_file, force_type="seq")
  seq_in.check_file_type("seq")
  sequences = seq_in.file_object

  cif_block = pdb_hierarchy.as_cif_block_with_sequence(
    sequences, crystal_symmetry=pdb_input.crystal_symmetry(),
    alignment_params=params)
  block_name = list(cif_model.keys())[0]

  def float_or_none(string):
    try: return float(string)
    except TypeError: return None

  d_min_from_cif = float_or_none(cif_block.get('_refine.ls_d_res_high'))
  d_max_from_cif = float_or_none(cif_block.get('_refine.ls_d_res_low'))
  # XXX maybe the values from the CIF (i.e. those actually used in the refinement)
  # should override the input params?
  if params.high_resolution is not None:
    params.high_resolution = d_min_from_cif
  if params.low_resolution is not None:
    params.low_resolution = d_max_from_cif

  if params.input.unmerged_data is not None:
    from mmtbx.command_line import cc_star
    result = cc_star.run(params=params.input, out=out)
    cif_block.update(result.as_cif_block())

  cif_model[block_name].update(cif_block)
  cif_model[block_name].sort(key=category_sort_function)
  print("Writing updated CIF file:", file=out)
  print("  " + params.output.cif_file, file=out)
  with open(params.output.cif_file, "wb") as f:
    print(cif_model, file=f)
  return

if __name__ == '__main__':
  run(sys.argv[1:])
