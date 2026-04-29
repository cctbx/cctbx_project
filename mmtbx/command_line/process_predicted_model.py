"""Replace B-value field and split model into compact domains"""
# LIBTBX_SET_DISPATCHER_NAME mmtbx.process_predicted_model
# LIBTBX_SET_DISPATCHER_NAME phenix.process_predicted_model
from __future__ import division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import process_predicted_model

if __name__ == '__main__':
  run_program(program_class=process_predicted_model.Program)

