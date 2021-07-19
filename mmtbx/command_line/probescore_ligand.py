from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME molprobity.probescore_ligand

import sys
from libtbx.utils import multi_out, show_total_time
from iotbx.cli_parser import CCTBXParser

from mmtbx.programs import probescore_ligand

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

  if len(cli_parser.namespace.unknown) > 0:
    # print("What is unknown: %s" % cli_parser.namespace.unknown)
    # print("Curr selection: '%s'" % wf.atom_selection_program.inselection)
    wf.atom_selection_program.inselection = cli_parser.namespace.unknown
    cli_parser.namespace.unknown = []
  cli_parser.working_phil = cli_parser.master_phil.format(python_object=wf)

def run(args):
  # create parser
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class = probescore_ligand.Program,
    custom_process_arguments = custom_args_proc,
    logger=logger)

  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)

  task = probescore_ligand.Program(
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
