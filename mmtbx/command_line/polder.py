# LIBTBX_SET_DISPATCHER_NAME phenix.polder
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import polder

if __name__ == '__main__':
  run_program(program_class=polder.Program)
