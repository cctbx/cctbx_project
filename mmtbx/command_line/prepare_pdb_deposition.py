from __future__ import division, print_function

import os, sys

import iotbx.pdb
import mmtbx.model

from iotbx.file_reader import any_file
from libtbx.data_manager import DataManager
from libtbx.utils import multi_out, Sorry
from mmtbx.programs import prepare_pdb_deposition

# =============================================================================
# roll into command-line parser
def show_usage(logger=None):
  if (logger is None):
    logger = sys.stdout
  print(prepare_pdb_deposition.description +
        '\n mmtbx.prepare_pdb_depostion <model file> <sequence file>\n',
        file=logger)

# =============================================================================
def run(args):

  logger = multi_out() #logging.getLogger('main')
  logger.register('stdout', sys.stdout)

  # replace command-line parser
  # ---------------------------------------------------------------------------
  if (len(args) == 0):
    show_usage(logger=logger)
    sys.exit()

  input_objects = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=prepare_pdb_deposition.master_params,
    pdb_file_def='input.model_file',
    seq_file_def='input.sequence_file'
  )

  # get program settings
  params = input_objects.work.extract()

  # get files (will be handled by task.validate)
  # already
  # if (params.input.model_file is None):
  #   raise Sorry('One model file is required.')
  # if (params.input.sequence_file is None):
  #   raise Sorry('One sequence file is required.')

  pdb_input = iotbx.pdb.input(params.input.model_file)
  model = mmtbx.model.manager(model_input=pdb_input, log=logger)

  sequence = any_file(params.input.sequence_file)
  sequence.check_file_type('seq')
  sequence = sequence.file_object

  # construct data manager
  data_manager = DataManager()
  data_manager.add_model(params.input.model_file, model)
  data_manager.add_sequence(params.input.sequence_file, sequence)

  # ---------------------------------------------------------------------------
  # start program
  task = prepare_pdb_deposition.Program(data_manager, params, logger=logger)

  # validate inputs
  task.validate()

  # run program
  task.run()

# =============================================================================
if __name__ == '__main__':
  run(sys.argv[1:])





# =============================================================================
# old code
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


def old_run(args, out=None):
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
  block_name = cif_model.keys()[0]

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
  print >> out, "Writing updated CIF file:"
  print >> out, "  " + params.output.cif_file
  with open(params.output.cif_file, "wb") as f:
    print >> f, cif_model
  return
