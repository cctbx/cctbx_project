from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.compare_rama

from mmtbx.programs import compare_rama
from iotbx.cli_parser import run_program

run_program(compare_rama.Program)
