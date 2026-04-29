# LIBTBX_SET_DISPATCHER_NAME mmtbx.debugging.riding
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import debug_riding

if __name__ == "__main__":
  run_program(debug_riding.Program)
