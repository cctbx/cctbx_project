from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.comparama

from mmtbx.programs import comparama
from iotbx.cli_parser import run_program

run_program(comparama.Program)
