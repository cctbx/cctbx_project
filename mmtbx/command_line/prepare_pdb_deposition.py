from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import prepare_pdb_deposition

if __name__ == '__main__':
  run_program(program_class=prepare_pdb_deposition.Program)
