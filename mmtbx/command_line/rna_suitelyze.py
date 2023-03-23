# LIBTBX_SET_DISPATCHER_NAME phenix.rna_suitelyze
# LIBTBX_SET_DISPATCHER_NAME molprobity.rna_suitelyze

from __future__ import absolute_import, division, print_function
import sys
from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time
from mmtbx.programs import rna_suitelyze

def run(args, out=sys.stdout, quiet=False):
  # create parser
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class=rna_suitelyze.Program,
    logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = rna_suitelyze.Program(
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

if (__name__ == "__main__"):
  run(sys.argv[1:])
