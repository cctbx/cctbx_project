# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_emdb
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import fetch_emdb

if __name__ == '__main__':
  run_program(program_class=fetch_emdb.Program)
