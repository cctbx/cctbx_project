# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.qmi
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import qmi

if __name__ == '__main__':
  run_program(program_class=qmi.Program)
