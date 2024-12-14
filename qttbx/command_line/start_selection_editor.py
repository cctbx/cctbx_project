# LIBTBX_SET_DISPATCHER_NAME phenix.start_selection_editor
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_selection_editor
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qttbx.programs.start_selection_editor import SelectionProgram

if __name__ == '__main__':
  run_program(program_class=SelectionProgram)
