# LIBTBX_SET_DISPATCHER_NAME phenix.start_molstar
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_molstar
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qttbx.programs import start_molstar

if __name__ == '__main__':
  run_program(program_class=start_molstar.Program)
