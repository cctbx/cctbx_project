from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME mmtbx.find_reference
# LIBTBX_SET_DISPATCHER_NAME phenix.find_reference

from mmtbx.programs import find_reference
from iotbx.cli_parser import run_program

result = run_program(find_reference.Program)
