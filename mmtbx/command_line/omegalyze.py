"""Counts of cis prolines and twisted prolines"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.omegalyze
# LIBTBX_SET_DISPATCHER_NAME molprobity.omegalyze
###
###import sys
###from mmtbx.validation import omegalyze
###
###if __name__ == "__main__":
###  omegalyze.run(sys.argv[1:])

import sys

from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time #, null_out

from iotbx.cli_parser import run_program
from mmtbx.programs import omegalyze

# =============================================================================
def old_run(args):

  # create parser
  #logger = multi_out() #logging.getLogger('main')
  #logger.register('stdout', sys.stdout)
  #logger = null_out()
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)
  #only omegalyze output is sent to stdout for backward compatibility with
  #  MolProbity website

  parser = CCTBXParser(
    program_class=omegalyze.Program,
    logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = omegalyze.Program(
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
  # old_run(sys.argv[1:])
  run_program(program_class=omegalyze.Program, hide_parsing_output=True)

