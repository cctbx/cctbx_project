"""
Provides a command-line utility for fetching PDB files and their associated
reflection data.
"""

from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.fetch_pdb

from mmtbx.programs import fetch
from iotbx.cli_parser import run_program

result = run_program(fetch.Program)
