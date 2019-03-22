# -*- coding: utf-8 -*-
from __future__ import division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.development.ntc_validation

import os, sys

from iotbx.cli_parser import run_program
from mmtbx.programs.ntc_validation import Program

# =============================================================================

if (__name__ == '__main__'):
  results = run_program(program_class=Program)
