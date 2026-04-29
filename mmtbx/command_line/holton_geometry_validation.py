"""Score a model with Holton geometry score"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.holton_geometry_validation
# LIBTBX_SET_DISPATCHER_NAME mmtbx.holton_geometry_validation
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from mmtbx.programs import holton_geometry_validation
from iotbx.cli_parser import run_program

if __name__ == '__main__':
  run_program(program_class=holton_geometry_validation.Program)

