"""Run ChimieraX"""
# LIBTBX_SET_DISPATCHER_NAME phenix.start_chimerax
# LIBTBX_SET_DISPATCHER_NAME qttbx.start_chimerax
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from qttbx.programs import start_chimerax

if __name__ == '__main__':
  run_program(program_class=start_chimerax.Program)

