from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.fetch_pdb

import sys
from mmtbx.programs import fetch
from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time

def custom_args_proc(cli_parser):
  wf = cli_parser.working_phil.extract()
  if len(cli_parser.namespace.unknown) > 0:
    # print("What is unknown: %s" % cli_parser.namespace.unknown)
    # print("Curr selection: '%s'" % wf.fetch.pdb_ids)
    wf.fetch.pdb_ids = cli_parser.namespace.unknown
    cli_parser.namespace.unknown = []
  cli_parser.working_phil = cli_parser.master_phil.format(python_object=wf)

def run(args):
  # create parser
  logger = multi_out()
  logger.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class = fetch.Program,
    custom_process_arguments = custom_args_proc,
    logger=logger)

  namespace = parser.parse_args(args)

  print('Starting job', file=logger)
  print('='*79, file=logger)

  task = fetch.Program(
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
