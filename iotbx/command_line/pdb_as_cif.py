"""Convert PDB formatted model to mmCIF"""
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_as_cif
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import pdb_as_cif

if __name__ == '__main__':
  run_program(program_class=pdb_as_cif.Program)
