from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import ribbons

if __name__ == "__main__":
  run_program(ribbons.Program)
