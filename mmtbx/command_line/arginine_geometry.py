from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs.arginine_geometry import Program

if (__name__ == '__main__'):
  results = run_program(program_class=Program)
