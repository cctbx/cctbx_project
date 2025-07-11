from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.ss_idealization

from mmtbx.programs import ss_idealization
from iotbx.cli_parser import run_program

if __name__ == "__main__":
  run_program(ss_idealization.Program)
