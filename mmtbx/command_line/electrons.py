# LIBTBX_SET_DISPATCHER_NAME mmtbx.electrons
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.ligands import electrons

run_program(electrons.Program)
