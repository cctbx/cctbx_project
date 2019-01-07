from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.comparama

from mmtbx.programs import comparama
from iotbx.cli_parser import run_program

run_program(comparama.Program)
