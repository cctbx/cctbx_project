from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.fix_cablam

from iotbx.cli_parser import run_program
from mmtbx.programs import fix_cablam

run_program(fix_cablam.Program)