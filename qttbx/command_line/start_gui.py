# LIBTBX_SET_DISPATCHER_NAME phenix.start_gui
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_gui
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qttbx.programs import start_gui

if __name__ == '__main__':
  run_program(program_class=start_gui.Program)
