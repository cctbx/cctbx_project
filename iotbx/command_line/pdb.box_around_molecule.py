"""Make a box around molecule and call it the unit cell"""
from __future__ import absolute_import, division, print_function
from iotbx.cli_parser import run_program
from iotbx.programs import box_around_molecule

if __name__ == '__main__':
  run_program(program_class=box_around_molecule.Program)
