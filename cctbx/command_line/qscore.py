# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.qscore
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from cctbx.programs import qscore

if __name__ == '__main__':
  run_program(program_class=qscore.Program)
