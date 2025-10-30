"""Compute map correlation coefficient given input model and reflection data"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.real_space_correlation

from iotbx.cli_parser import run_program
from mmtbx.programs import real_space_correlation

if __name__ == '__main__':
  run_program(program_class=real_space_correlation.Program)

