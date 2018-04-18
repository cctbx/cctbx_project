from __future__ import division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.fix_cablam

import sys

from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time
from mmtbx.programs import fix_cablam


# =============================================================================
def run(args):

  # create parser
  logger = multi_out() #logging.getLogger('main')
  logger.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class=fix_cablam.Program,
    logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = fix_cablam.Program(
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
