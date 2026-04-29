from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

if __name__ == "__main__":
  run_program(reduce2.Program)
