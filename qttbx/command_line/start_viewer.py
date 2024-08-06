# LIBTBX_SET_DISPATCHER_NAME phenix.start_viewer
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_viewer
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qttbx.programs import start_viewer

if __name__ == '__main__':
  run_program(program_class=start_viewer.Program)
