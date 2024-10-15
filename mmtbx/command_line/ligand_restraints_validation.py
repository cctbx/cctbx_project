# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.ligand_restraints_validation
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs.ligand_restraints_validation import Program

if (__name__ == '__main__'):
  results = run_program(program_class=Program)
