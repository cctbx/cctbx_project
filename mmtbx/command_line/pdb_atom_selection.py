from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_atom_selection

import sys
from libtbx.utils import multi_out, show_total_time
from iotbx.cli_parser import CCTBXParser

from mmtbx.programs import atom_selection

def custom_args_proc(cli_parser):
  """
  This is going to put inselection positional argument into phil scope of the
  program.
  Useful things are:
  cli_parser.namespace
  cli_parser.data_manager
  cli_parser.namespace.unknown - where these argument is going to be
  cli_parser.working_phil
  """
  wf = cli_parser.working_phil.extract()

  # Since we have phil parameters for these as well, we want to make
  # sure that we don't overwrite them if nothing was provided via
  # command-line.
  if cli_parser.namespace.cryst1_replacement_buffer_layer is not None:
    wf.atom_selection_program.cryst1_replacement_buffer_layer = \
        cli_parser.namespace.cryst1_replacement_buffer_layer
  if cli_parser.namespace.write_pdb_file is not None:
    wf.atom_selection_program.write_pdb_file = \
        cli_parser.namespace.write_pdb_file

  if len(cli_parser.namespace.unknown) > 0:
    # print("What is unknown: %s" % cli_parser.namespace.unknown)
    # print("Curr selection: '%s'" % wf.atom_selection_program.inselection)
    wf.atom_selection_program.inselection = cli_parser.namespace.unknown
    cli_parser.namespace.unknown = []
  cli_parser.working_phil = cli_parser.master_phil.format(python_object=wf)


def run(args):
  # create parser
  logger = multi_out() #logging.getLogger('main')
  logger.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class = atom_selection.Program,
    custom_process_arguments = custom_args_proc,
    logger=logger)
  #
  # Stick in additional keyword args
  # They going to end up in namespace.cryst1_replacement_buffer_layer etc
  # file_name is obsolet, parser takes care of it putting it in data manager
  # inselections is going to be handled by custom_args_proc function
  # because it is intended to be positional argument
  #
  # !!! This is done here for legacy support and illustrative purposes.
  # Don't do it anywhere else, since we strive to have the same
  # command-line flags across all programs, like --overwrite etc.
  parser.add_argument(
      "--cryst1-replacement-buffer-layer",
      action="store",
      type=float,
      help="replace CRYST1 with pseudo unit cell covering the selected"
        " atoms plus a surrounding buffer layer",
      default=None)
  parser.add_argument(
      "--write-pdb-file", "--write_pdb_file", "--write_pdb-file", "--write-pdb_file",
      action="store",
      help="write selected atoms to new PDB file",
      default=None)

  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)

  task = atom_selection.Program(
      parser.data_manager, parser.working_phil.extract(),
      master_phil=parser.master_phil,
      logger=logger)

  # validate inputs
  task.validate()

  # run program
  task.run()

  # stop timer
  print('', file=logger)
  print('='*79, file=logger)
  print('Job complete', file=logger)
  show_total_time(out=logger)

if (__name__ == "__main__"):
  run(sys.argv[1:])
