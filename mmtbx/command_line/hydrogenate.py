# LIBTBX_SET_DISPATCHER_NAME mmtbx.debugging.hydrogenate
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import hydrogenate

if __name__ == "__main__":
  run_program(hydrogenate.Program)
