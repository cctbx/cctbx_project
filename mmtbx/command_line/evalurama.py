from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.evalurama

from mmtbx.programs import evalurama
from iotbx.cli_parser import run_program

run_program(evalurama.Program)
