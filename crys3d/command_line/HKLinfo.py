"""Summarize reflection file"""
# LIBTBX_SET_DISPATCHER_NAME phenix.HKLinfo
# LIBTBX_SET_DISPATCHER_NAME cctbx.HKLinfo
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from crys3d.programs import HKLinfo

if __name__ == '__main__':
  run_program(program_class=HKLinfo.Program)
