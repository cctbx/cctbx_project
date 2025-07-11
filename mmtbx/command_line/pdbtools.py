"""Manipulate PDB files"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.pdbtools

from iotbx.cli_parser import run_program
from mmtbx.programs import pdbtools

run_program(pdbtools.Program)

