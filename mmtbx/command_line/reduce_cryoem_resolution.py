# LIBTBX_SET_DISPATCHER_NAME phenix.reduce_cryoem_resolution
# LIBTBX_SET_DISPATCHER_NAME mmtbx.reduce_cryoem_resolution
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import reduce_cryoem_resolution

if __name__ == '__main__':
  run_program(program_class=reduce_cryoem_resolution.Program)
