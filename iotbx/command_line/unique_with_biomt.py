from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME iotbx.unique_with_biomt
# LIBTBX_SET_DISPATCHER_NAME phenix.unique_with_biomt

from iotbx.programs import unique_with_biomt
from iotbx.cli_parser import run_program

result = run_program(unique_with_biomt.Program)
