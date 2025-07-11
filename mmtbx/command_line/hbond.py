"""Tool to find H bonds in an atomic model"""
from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.hbond

from mmtbx.programs import hbond
from iotbx.cli_parser import run_program

if __name__ == "__main__":
  run_program(hbond.Program)

