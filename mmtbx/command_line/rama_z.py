"""Tool to calculate Rama-Z score. Validation of Ramachandran plot"""
from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME mmtbx.rama_z
# LIBTBX_SET_DISPATCHER_NAME phenix.rama_z

from mmtbx.programs import rama_z
from iotbx.cli_parser import run_program

if __name__ == "__main__":
  result = run_program(rama_z.Program)

