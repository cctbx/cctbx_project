from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.electrons

from iotbx.cli_parser import run_program
from mmtbx.ligands import electrons

if __name__ == "__main__":
  run_program(electrons.Program)
