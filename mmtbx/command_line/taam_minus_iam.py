from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.TAAM_minus_IAM
# LIBTBX_SET_DISPATCHER_NAME mmtbx.TAAM_minus_IAM
from iotbx.cli_parser import run_program
from mmtbx.programs import taam_minus_iam

if __name__ == '__main__':
  run_program(program_class=taam_minus_iam.Program)
