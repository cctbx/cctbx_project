from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.secondary_structure_validation

from mmtbx.programs import ss_validation
from iotbx.cli_parser import run_program

run_program(ss_validation.Program)
