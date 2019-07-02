from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.hbond

from mmtbx.programs import hbond
from iotbx.cli_parser import run_program

run_program(hbond.Program)
