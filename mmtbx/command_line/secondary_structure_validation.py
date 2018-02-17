from __future__ import division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.secondary_structure_validation

import sys

from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time
from mmtbx.programs import ss_validation


# =============================================================================
def run(args):

  # create parser
  logger = multi_out() #logging.getLogger('main')
  logger.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class=ss_validation.Program,
    logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = ss_validation.Program(
    parser.data_manager, parser.working_phil.extract(), logger=logger)

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
  run(sys.argv[1:])
