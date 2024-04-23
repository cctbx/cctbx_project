from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.water_b_factors

from iotbx.cli_parser import run_program
from mmtbx.programs.water_b_factors import Program

if (__name__ == '__main__'):
  results = run_program(program_class=Program)
