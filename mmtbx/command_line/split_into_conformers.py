"""Split a model with alternative conformations into complete conformers"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.split_into_conformers

from mmtbx.programs import split_into_conformers
from iotbx.cli_parser import run_program

if __name__ == "__main__":
  run_program(split_into_conformers.Program)
