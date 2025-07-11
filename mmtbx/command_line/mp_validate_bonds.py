"""Validate bond lengths"""
# LIBTBX_SET_DISPATCHER_NAME phenix.mp_validate_bonds
# LIBTBX_SET_DISPATCHER_NAME molprobity.mp_validate_bonds

from __future__ import absolute_import, division, print_function
import sys
from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time
from mmtbx.programs import mp_validate_bonds
from iotbx.cli_parser import run_program

def old_run(args, out=sys.stdout, quiet=False):
  # create parser
  logger = multi_out()
  logger.register('stderr', sys.stderr)
  logger2 = multi_out()
  logger2.register('stdout', sys.stdout)

  parser = CCTBXParser(
    program_class=mp_validate_bonds.Program,
    logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = mp_validate_bonds.Program(
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
  #run(sys.argv[1:])
  run_program(program_class=mp_validate_bonds.Program, hide_parsing_output=True)

