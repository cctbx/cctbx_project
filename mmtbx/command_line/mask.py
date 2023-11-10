from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.mask

from iotbx.cli_parser import run_program
from mmtbx.programs import mask

run_program(mask.Program)

