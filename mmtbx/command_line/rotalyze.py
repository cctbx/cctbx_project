"""Validate protein sidechain rotamers"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.rotalyze
# LIBTBX_SET_DISPATCHER_NAME molprobity.rotalyze
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

import sys

from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time
from mmtbx.programs import rotalyze
from iotbx.cli_parser import run_program

# =============================================================================
def old_run(args):

  # create parser
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class=rotalyze.Program,
    logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = rotalyze.Program(
    parser.data_manager, parser.working_phil.extract(), logger=logger2)

  # validate inputs
  task.validate()

  # run program
  task.run()

  # stop timer
  print('', file=logger)
  print('='*79, file=logger)
  print('Job complete', file=logger)
  show_total_time(out=logger)

# =============================================================================
if __name__ == '__main__':
  #run(sys.argv[1:])
  run_program(program_class=rotalyze.Program, hide_parsing_output=True)

