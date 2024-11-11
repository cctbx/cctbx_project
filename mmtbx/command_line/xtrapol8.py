# LIBTBX_SET_DISPATCHER_NAME phenix.development.xtrapol8
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import xtrapol8

if __name__ == '__main__':
  run_program(program_class=xtrapol8.Program)
