"""Identify unusual AlphaFold conformations"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.barbed_wire_analysis
# LIBTBX_SET_DISPATCHER_NAME molprobity.barbed_wire_analysis

from iotbx.cli_parser import run_program
from mmtbx.programs import barbed_wire_analysis

if __name__ == "__main__":
  run_program(barbed_wire_analysis.Program)

