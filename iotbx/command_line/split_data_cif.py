"""Split an mmCIF file"""
from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME iotbx.split_data_cif

from iotbx.programs import split_data_cif
from iotbx.cli_parser import run_program

result = run_program(split_data_cif.Program)
